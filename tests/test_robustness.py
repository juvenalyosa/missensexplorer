"""Failure handling, provenance, local VEP plumbing and streaming input (all offline)."""

import gzip
import importlib
import json
import stat
import sys
from pathlib import Path

import httpx
import pytest

from missensexplorer.api import ApiError, Http, OfflineError, SourceStatus
from missensexplorer.inputs import read_file
from missensexplorer.report import run_info, to_bytes
from missensexplorer.vep_local import LocalVep, to_default_format

P = importlib.import_module("missensexplorer.predict")  # the package re-exports a function with the same name
DATA = Path(__file__).parent / "data"
BRAF = json.loads((DATA / "vep_braf_v600e_region.json").read_text())


# --------------------------------------------------------------------------- HTTP


def test_offline_host_fails_fast(monkeypatch):
    monkeypatch.setattr(Http, "_sleep", staticmethod(lambda attempt: None))
    calls = []

    def handler(request):
        calls.append(request.url.host)
        raise httpx.ConnectError("no route to host", request=request)

    http = Http()
    http.client = httpx.Client(transport=httpx.MockTransport(handler))
    with pytest.raises(OfflineError):
        http.request("GET", "https://rest.ensembl.org/info/ping")
    n = len(calls)
    assert n == http.transport_retries
    with pytest.raises(OfflineError):  # host is now marked down: no new attempts
        http.request("GET", "https://rest.ensembl.org/info/software")
    assert len(calls) == n


def test_host_that_keeps_failing_is_skipped(monkeypatch):
    pauses = []
    monkeypatch.setattr(Http, "_pause", staticmethod(pauses.append))
    calls = []

    def handler(request):
        calls.append(request.url.path)
        if request.url.path == "/ok":
            return httpx.Response(200, json={})
        return httpx.Response(429 if "files" in request.url.path else 503)

    http = Http()
    http.client = httpx.Client(transport=httpx.MockTransport(handler))
    for i in range(2):
        with pytest.raises(ApiError, match="rate limited"):
            http.request("GET", f"https://alphafold.ebi.ac.uk/files/{i}.json")
    assert http.request("GET", "https://alphafold.ebi.ac.uk/ok").status_code == 200  # success resets the count
    for i in range(3):
        with pytest.raises(ApiError):
            http.request("GET", f"https://alphafold.ebi.ac.uk/files/x{i}.json")
    n = len(calls)
    assert len(pauses) == 5 * (http.max_retries - 1)  # no pause after the last attempt
    with pytest.raises(OfflineError, match="keeps failing: HTTP 429"):
        http.request("GET", "https://alphafold.ebi.ac.uk/files/later.json")
    assert len(calls) == n  # no new attempt once the host is marked down
    with pytest.raises(ApiError, match="HTTP 503"):  # other hosts are unaffected
        http.request("GET", "https://rest.uniprot.org/x")
    assert n == 5 * http.max_retries + 1  # each failed request tried max_retries times


def test_source_status_summary():
    st = SourceStatus()
    st.use("UniProt", 3)
    st.use("ClinVar", 2)
    st.fail("ClinVar", "HTTP 503", ["VCV1"])
    st.skip("AlphaFold DB", "offline")
    s = st.summary()
    assert s["UniProt"] == "ok"
    assert s["ClinVar"].startswith("failed (1 items)")
    assert s["AlphaFold DB"] == "skipped: offline"
    assert st.failed("ClinVar", "VCV1")
    assert not st.failed("ClinVar", "VCV2")
    assert not st.failed("UniProt")


# --------------------------------------------------------------------------- pipeline with fake services


class FakeVep:
    """Base annotation always works; the predictor plugins fail (as when a plugin times out)."""

    def __init__(self, http, cache, assembly, workers=4, release=None, status=None):
        self.status = status

    def annotate(self, queries, params, progress=None, stage="", source="Ensembl VEP"):
        out = {}
        for kind, qs in queries.items():
            for q in set(qs):
                if "AlphaMissense" in params:
                    out[(kind, q)] = {"error": "VEP request failed: timeout", "service_error": True}
                    self.status.fail(source, "timeout", [q])
                else:
                    out[(kind, q)] = {"results": BRAF}
        return out


@pytest.fixture
def fake_services(monkeypatch):
    monkeypatch.setattr(P, "ensembl_release", lambda http, cache, asm: 116)
    monkeypatch.setattr(P, "_match_rate", lambda *a, **k: (1.0, {}))
    monkeypatch.setattr(P, "Vep", FakeVep)

    def clinvar_down(http, cache, ids, status=None):
        status.use("ClinVar", len(ids))
        status.fail("ClinVar", "HTTP 503", ids)
        return {}

    def uniprot_down(http, cache, accs, workers=4, status=None):
        status.fail("UniProt", "HTTP 500", accs)
        return {}

    def clinvar_search_down(http, cache, windows, assembly, status=None):
        status.fail("ClinVar", "HTTP 503", [f"{assembly}:{c}:{a}-{b}" for c, a, b in windows])
        return {}

    def gnomad_down(http, cache, variants, assembly, status=None, progress=None):
        status.fail("gnomAD", "HTTP 503", [f"{c}-{p}-{r}-{a}" for c, p, r, a in variants])
        return {}

    monkeypatch.setattr(P.cv, "fetch", clinvar_down)
    monkeypatch.setattr(P.cv, "search_windows", clinvar_search_down)
    monkeypatch.setattr(P.gnomad, "lookup", gnomad_down)
    monkeypatch.setattr(P.prot, "fetch_uniprot", uniprot_down)
    monkeypatch.setattr(P.prot, "fetch_alphafold", lambda *a, **k: {})
    monkeypatch.setattr(P.prot, "fetch_protein_sequences", lambda *a, **k: {})


def test_predictor_failure_keeps_annotation(fake_services):
    res = P.predict(text="7 140753336 A T", options=P.PredictOptions(cache=False))
    assert len(res.table) == 1
    assert len(res.errors) == 0
    r = res.table.iloc[0]
    assert r["protein_change"] == "V600E"
    assert "predictors unavailable" in r["notes"]
    assert "ClinVar unavailable" in r["notes"]
    assert "UniProt unavailable" in r["notes"]
    assert res.sources["Ensembl VEP + predictors"].startswith("failed")
    assert res.sources["ClinVar"].startswith("failed")
    assert any("ClinVar did not respond" in w for w in res.warnings)


def test_two_pass_predictor_failure_keeps_annotation(fake_services):
    # > 200 unique variants -> fast pass + predictor pass; the predictor pass fails
    text = "\n".join(f"7 {140753336 + i} A T" for i in range(205))
    res = P.predict(text=text, options=P.PredictOptions(cache=False, clinvar=False, protein_context=False))
    assert len(res.errors) == 0
    assert len(res.table) == 205
    assert res.table["notes"].str.contains("predictors unavailable").all()
    assert res.sources["VEP predictor plugins"].startswith("failed")


def test_run_record_and_exports(fake_services, tmp_path):
    vcf = tmp_path / "s.vcf"
    vcf.write_text(
        "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n7\t140753336\t.\tA\tT\t.\tPASS\t.\n"
    )
    res = P.predict(files=[str(vcf)], options=P.PredictOptions(cache=False))
    run = res.run
    assert run["ensembl_release"] == 116
    assert "release 116" in run["annotation"]
    assert run["inputs"]["files"][0]["name"] == "s.vcf"
    assert len(run["inputs"]["files"][0]["sha256"]) == 64
    assert run["options"]["vep"] == "rest"
    info = run_info(res)
    assert "source: ClinVar" in set(info["item"])
    assert json.loads(to_bytes(res, "json"))["run"]["tool"].startswith("MissensExplorer")
    assert len(to_bytes(res, "xlsx")) > 1000


def test_offline_requires_local_vep():
    with pytest.raises(ValueError, match="offline mode needs the local VEP"):
        P.predict(text="7 140753336 A T", options=P.PredictOptions(offline=True))


# --------------------------------------------------------------------------- local VEP


def test_to_default_format():
    assert to_default_format("7 140753336 . A T . . .", "x") == "7\t140753336\t140753336\tA/T\t+\tx"
    assert to_default_format("7 117559590 . ATCT A . . .", "x") == "7\t117559591\t117559593\tTCT/-\t+\tx"
    assert to_default_format("7 100 . A AT . . .", "x") == "7\t101\t100\t-/T\t+\tx"
    assert to_default_format("7 100 100 A/- 1", "x") == "7\t100\t100\tA/-\t1\tx"


@pytest.fixture
def fake_vep_cmd(tmp_path):
    """A stand-in for `vep`: echoes the BRAF fixture for each stdin line, with its identifier."""
    script = tmp_path / "fakevep.py"
    script.write_text(
        "import json, sys\n"
        f"fixture = json.load(open({str(DATA / 'vep_braf_v600e_region.json')!r}))[0]\n"
        "if '--help' in sys.argv:\n    print('ensembl-vep : 115.0'); sys.exit(0)\n"
        "for line in sys.stdin:\n"
        "    f = line.rstrip('\\n').split('\\t')\n"
        "    if len(f) < 6: continue\n"
        "    if f[0] == 'BAD': continue\n"
        "    item = dict(fixture); item['id'] = f[5]; item['input'] = line.strip()\n"
        "    print(json.dumps(item))\n"
    )
    script.chmod(script.stat().st_mode | stat.S_IEXEC)
    return f"{sys.executable} {script}"


def test_local_vep_backend(fake_vep_cmd, monkeypatch):
    monkeypatch.setattr(P, "ensembl_release", lambda http, cache, asm: 116)
    monkeypatch.setattr(P, "_match_rate", lambda *a, **k: (1.0, {}))
    opts = P.PredictOptions(vep="local", vep_command=fake_vep_cmd, offline=True, cache=False)
    res = P.predict(text="7 140753336 A T\nBRAF V600E", options=opts)
    assert list(res.table["protein_change"]) == ["V600E"]
    assert "HGVS and rsIDs need the online mode" in res.errors.iloc[0]["error"]
    assert res.run["annotation"] == "local VEP (ensembl-vep 115.0)"
    assert res.sources["ClinVar"].startswith("skipped")


def test_local_vep_missing_command_is_reported():
    st = SourceStatus()
    lv = LocalVep(P.NullCache(), command="/nonexistent/vep", status=st)
    out = lv.annotate({"region": ["7 140753336 . A T . . ."]})
    val = out[("region", "7 140753336 . A T . . .")]
    assert val["service_error"]
    assert "could not run local VEP" in val["error"]
    assert st.failed("Ensembl VEP (local)")


# --------------------------------------------------------------------------- streaming input


def test_vcf_gz_is_streamed(tmp_path, monkeypatch):
    p = tmp_path / "big.vcf.gz"
    with gzip.open(p, "wt") as fh:
        fh.write("##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n")
        for i in range(2000):
            fh.write(f"1\t{1000 + i}\t.\tA\tG\t.\tPASS\t.\tGT\t0/1\n")
    monkeypatch.setattr(Path, "read_bytes", lambda self: (_ for _ in ()).throw(AssertionError("loaded whole file")))
    parsed = read_file(str(p))
    assert len(parsed.variants) == 2000
    assert parsed.variants[0].zygosity == "heterozygous"


def test_vep_retries_truncated_response():
    from missensexplorer.vep import Vep

    calls = []

    class FakeHttp:
        def request(self, method, url, **kw):
            calls.append(kw["json"]["variants"])
            req = httpx.Request("POST", url)
            if len(calls) == 1:  # first answer is cut off mid-JSON
                return httpx.Response(200, text=json.dumps(BRAF)[:500], request=req)
            return httpx.Response(200, text=json.dumps(BRAF), request=req)

    vep = Vep(FakeHttp(), P.NullCache(), "GRCh38", workers=1)
    q = BRAF[0]["input"]
    out = vep.annotate({"region": [q]}, {"hgvs": 1})
    assert "results" in out[("region", q)]
    assert len(calls) == 2


def test_votes_accept_table_rows_with_nan():
    nan = float("nan")
    row = {
        "alphamissense": 0.99,
        "alphamissense_class": "likely_pathogenic",
        "revel": nan,
        "cadd_phred": nan,
        "eve_class": nan,
        "sift": nan,
        "polyphen": None,
    }
    dmg, ben, _ = P._votes(row)
    assert dmg == ["AlphaMissense"]
    assert ben == []


def test_one_unreadable_variant_does_not_sink_its_batch():
    from missensexplorer.vep import Vep

    good = BRAF[0]["input"]
    bad = "1 12345 . A G . . ."

    class FakeHttp:
        def request(self, method, url, **kw):
            req = httpx.Request("POST", url)
            if bad in kw["json"]["variants"]:
                return httpx.Response(200, text='[{"input": "1 12345 . A G . . .", "x": ', request=req)
            return httpx.Response(200, text=json.dumps(BRAF), request=req)

    st = SourceStatus()
    vep = Vep(FakeHttp(), P.NullCache(), "GRCh38", workers=1, status=st)
    out = vep.annotate({"region": [good, bad]}, {"hgvs": 1})
    assert "results" in out[("region", good)]
    assert out[("region", bad)]["service_error"]
    assert st.failed("Ensembl VEP", bad)


def test_an_unreachable_ensembl_with_a_filled_cache_uses_the_cache(monkeypatch, tmp_path):
    from types import SimpleNamespace

    from missensexplorer.api import OfflineError
    from missensexplorer.cache import Cache

    def offline(*a, **k):
        raise OfflineError("no network")

    monkeypatch.setattr(P, "ensembl_release", offline)
    cache = Cache(tmp_path / "c.sqlite")
    ctx = SimpleNamespace(net=True, local=False, cache=cache, http=None, warnings=[])
    with pytest.raises(ConnectionError):  # nothing cached: no answers to use
        P._check_connection(ctx)
    cache.set("ensembl_release|GRCh38", 116)
    P._check_connection(ctx)  # a rerun from a filled cache goes on, and says so
    assert ctx.warnings == ["Ensembl cannot be reached: only answers already in the cache are used"]
    cache.close()


def test_a_rerun_from_the_cache_without_the_internet_gives_the_same_rows(monkeypatch, tmp_path):
    """A run fills the cache; two days later, with Ensembl unreachable, the same input gives the same rows."""
    import time as _time

    import missensexplorer.cache as cache_mod
    from missensexplorer.api import Http

    braf = json.loads((DATA / "vep_braf_v600e_region.json").read_text())
    online = {"on": True}

    def handler(request):
        if not online["on"]:
            raise httpx.ConnectError("no network", request=request)
        if request.url.path.endswith("/info/software"):
            return httpx.Response(200, json={"release": 116})
        if "/vep/" in request.url.path:
            return httpx.Response(200, json=braf)
        return httpx.Response(404, json={})

    class MockHttp(Http):
        def __init__(self, *a, **k):
            super().__init__(*a, max_retries=1, transport_retries=1, max_wait=0, **k)
            self.client = httpx.Client(transport=httpx.MockTransport(handler))

    monkeypatch.setattr(P, "Http", MockHttp)
    opts = {"cache_path": str(tmp_path / "c.sqlite"), "clinvar": False, "protein_context": False,
            "sequences": False, "gnomad_check": False}  # fmt: skip
    first = P.predict(text="7 140753336 A T", options=P.PredictOptions(**opts))
    assert len(first.table) == 1
    online["on"] = False
    later = _time.time() + 2 * 86400  # the cached Ensembl release (kept a day) has expired
    monkeypatch.setattr(cache_mod.time, "time", lambda: later)
    again = P.predict(text="7 140753336 A T", options=P.PredictOptions(**opts))
    assert "Ensembl cannot be reached: only answers already in the cache are used" in again.warnings
    assert again.run["cache_only"]
    assert again.run["ensembl_release"] == 116
    cols = ["gene", "protein_change", "prediction", "acmg_pp3_bp4"]
    assert again.table[cols].equals(first.table[cols])
