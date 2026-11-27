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
    monkeypatch.setattr(P, "ensembl_release", lambda http, cache, asm, **k: 116)
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
    monkeypatch.setattr(P, "ensembl_release", lambda http, cache, asm, **k: 116)
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


CACHE_ONLY = (
    "Ensembl cannot be reached: only Ensembl answers already in the cache are used, whatever their age; inputs "
    "without them fail"
)


class _Ping:
    """A stand-in client for the second request that tells whether Ensembl is down."""

    def __init__(self, error: ApiError | None) -> None:
        self.error, self.down, self.asked = error, {}, 0

    def get_json(self, url, **kw):
        self.asked += 1
        if self.error:
            raise self.error
        return {"ping": 1}


@pytest.mark.parametrize(
    ("error", "ping", "frozen"),
    [
        (OfflineError("no network"), None, True),  # no connection
        (
            ApiError("HTTP 503 after the retries"),
            ApiError("HTTP 503"),
            True,
        ),  # an outage that a second request confirms
        (ApiError("read timeout after the retries"), OfflineError("down"), True),
        (ApiError("HTTP 503 after the retries"), None, False),  # a passing error: the server answers again
        (ApiError("HTTP 403", 403), None, False),  # a client error: the server is up
        (ApiError("HTTP 404", 404), ApiError("HTTP 503"), False),
    ],
)
def test_the_cache_alone_is_used_only_when_ensembl_is_down(monkeypatch, tmp_path, error, ping, frozen):
    from types import SimpleNamespace

    from missensexplorer.cache import Cache

    def failing(*a, **k):
        raise error

    monkeypatch.setattr(P, "ensembl_release", failing)
    cache = Cache(tmp_path / "c.sqlite")
    cache.set("ensembl_release|GRCh38", 116)
    http = _Ping(ping)
    ctx = SimpleNamespace(
        net=True, local=False, cache=cache, http=http, assembly="GRCh38", warnings=[], cache_only=False
    )
    P._check_connection(ctx)
    assert ctx.cache_only is frozen
    assert cache.frozen is frozen
    assert ctx.warnings == ([CACHE_ONLY] if frozen else [])
    assert ("rest.ensembl.org" in http.down) is frozen  # in a run from the cache, Ensembl is not asked again
    assert http.asked == (0 if isinstance(error, OfflineError) or error.status else 1)
    cache.close()


def test_an_unreachable_ensembl_without_a_cache_for_the_build(monkeypatch, tmp_path):
    from types import SimpleNamespace

    from missensexplorer.cache import Cache

    def failing(*a, **k):
        raise failing.error

    monkeypatch.setattr(P, "ensembl_release", failing)
    cache = Cache(tmp_path / "c.sqlite")
    cache.set("ensembl_release|GRCh37", 75)  # a cache of the other build holds no answer for this one
    ctx = SimpleNamespace(net=True, local=False, cache=cache, http=_Ping(ApiError("HTTP 503")), assembly="GRCh38",
                          warnings=[], cache_only=False)  # fmt: skip
    failing.error = OfflineError("no network")
    with pytest.raises(ConnectionError):  # nothing to use
        P._check_connection(ctx)
    failing.error = ApiError("HTTP 503 after the retries")
    P._check_connection(ctx)  # a server outage without a cache: the run goes on and its requests fail one by one
    assert not ctx.cache_only
    assert not cache.frozen
    assert ctx.warnings == []
    cache.close()


def test_a_live_answer_is_never_keyed_to_a_cached_release(monkeypatch, tmp_path):
    from types import SimpleNamespace

    from missensexplorer.cache import Cache

    def failing(*a, **k):
        raise ApiError("HTTP 503 after the retries")

    monkeypatch.setattr(P, "ensembl_release", failing)
    cache = Cache(tmp_path / "c.sqlite")
    cache.set("ensembl_release|GRCh38", 116)
    ctx = SimpleNamespace(net=True, cache=cache, assembly="GRCh38", cache_only=False, http=None)
    assert P._release(ctx) is None  # the server's release is unknown: its answers are not filed under 116
    ctx.cache_only = True
    assert P._release(ctx) == 116  # a run from the cache alone reads the answers filed under 116
    cache.close()


def _ensembl_mock(monkeypatch, tmp_path):
    """A mocked Ensembl whose behaviour the test sets per endpoint; returns the state, the options and the log."""
    import missensexplorer.cache as cache_mod
    from missensexplorer.api import Http

    braf = json.loads((DATA / "vep_braf_v600e_region.json").read_text())
    state = {"software": "ok", "ping": "ok", "vep": "ok", "release": 116, "asked": []}

    def answer(kind, request):
        if kind == "none":
            raise httpx.ConnectError("no network", request=request)
        if kind == "timeout":
            raise httpx.ReadTimeout("no answer", request=request)
        if isinstance(kind, int):
            return httpx.Response(kind, json={})
        return None

    def handler(request):
        path = request.url.path
        state["asked"].append(path)
        if path.endswith("/info/software"):
            return answer(state["software"], request) or httpx.Response(200, json={"release": state["release"]})
        if path.endswith("/info/ping"):
            return answer(state["ping"], request) or httpx.Response(200, json={"ping": 1})
        if "/vep/" in path:
            return answer(state["vep"], request) or httpx.Response(200, json=braf)
        return httpx.Response(404, json={})

    class MockHttp(Http):
        def __init__(self, *a, **k):
            super().__init__(*a, max_retries=1, transport_retries=1, max_wait=0, **k)
            self.client = httpx.Client(transport=httpx.MockTransport(handler))

    monkeypatch.setattr(P, "Http", MockHttp)
    clock = {"t": None}
    real = cache_mod.time.time
    monkeypatch.setattr(cache_mod.time, "time", lambda: clock["t"] or real())
    opts = {"cache_path": str(tmp_path / "c.sqlite"), "clinvar": False, "protein_context": False,
            "sequences": False, "gnomad_check": False}  # fmt: skip
    return state, opts, clock


def _keys(path, prefix):
    import sqlite3

    with sqlite3.connect(path) as db:
        return {k for (k,) in db.execute("SELECT k FROM kv WHERE k LIKE ?", (prefix + "%",))}


def test_a_rerun_from_the_cache_without_the_internet_gives_the_same_rows(monkeypatch, tmp_path):
    """A run fills the cache; a month later, with Ensembl down, the same input gives the same rows."""
    import time as _time

    state, opts, clock = _ensembl_mock(monkeypatch, tmp_path)
    first = P.predict(text="7 140753336 A T", options=P.PredictOptions(**opts))
    assert len(first.table) == 1
    # the mock has no genome sequence: the REF allele cannot be checked, and the run says so
    assert "REF allele not checked against the genome for 1 input(s) (Ensembl sequence unavailable)" in first.warnings
    clock["t"] = _time.time() + 31 * 86400  # past the cache's 30-day limit: only a frozen cache still answers
    cols = ["gene", "protein_change", "prediction", "acmg_pp3_bp4"]
    for down in ("none", 503, "timeout", 429):  # no connection; an outage; no answer; rate limits, after the retries
        state.update(software=down, ping=down, vep=down, asked=[])
        again = P.predict(text="7 140753336 A T", options=P.PredictOptions(**opts))
        assert CACHE_ONLY in again.warnings
        assert again.run["cache_only"]
        assert again.run["ensembl_release"] == 116
        assert again.table[cols].equals(first.table[cols]), down
        assert not any("/vep/" in a for a in state["asked"])  # every Ensembl answer came from the cache
        # an input whose answers are not cached fails, and Ensembl is not asked for it
        new = P.predict(text="7 140753336 A C", options=P.PredictOptions(**opts))
        assert new.table.empty
        assert len(new.errors) == 1
        assert not any("/vep/" in a for a in state["asked"])


def test_a_passing_ensembl_error_does_not_freeze_the_run(monkeypatch, tmp_path):
    """One failed release lookup on a working Ensembl: the run stays live and files nothing under the old release."""
    import time as _time

    state, opts, clock = _ensembl_mock(monkeypatch, tmp_path)
    P.predict(text="7 140753336 A T", options=P.PredictOptions(**opts))
    old = _keys(opts["cache_path"], "vep|e116|")
    assert old
    clock["t"] = _time.time() + 2 * 86400
    state["release"] = 117  # Ensembl has moved on; its release lookup fails, the rest answers
    for software in (503, "timeout", 403, 404):
        state.update(software=software, asked=[])
        run = P.predict(text="7 140753336 A T", options=P.PredictOptions(predictors=False, **opts))
        assert not run.run["cache_only"], software
        assert CACHE_ONLY not in run.warnings
        assert run.run["ensembl_release"] is None  # unknown, never the cached 116
        if software == 503:  # asked live (later runs find the answer filed under the unknown release)
            assert any("/vep/" in a for a in state["asked"])
        assert _keys(opts["cache_path"], "vep|e116|") == old  # no live answer filed under the old release
    assert _keys(opts["cache_path"], "vep|e?|")


def test_a_truncation_across_an_exon_intron_boundary_has_uncertain_nmd(monkeypatch):
    import copy

    braf = copy.deepcopy(BRAF)
    for tc in braf[0]["transcript_consequences"]:
        tc.update(consequence_terms=["frameshift_variant"], exon="15/18", intron="15/17")

    class Vep:
        def __init__(self, *a, **k):
            pass

        def annotate(self, queries, params, progress=None, stage="", source="Ensembl VEP"):
            return {(kind, q): {"results": braf} for kind, qs in queries.items() for q in set(qs)}

    monkeypatch.setattr(P, "ensembl_release", lambda http, cache, asm, **k: 116)
    monkeypatch.setattr(P, "_match_rate", lambda *a, **k: (1.0, {}))
    monkeypatch.setattr(P, "Vep", Vep)
    opts = P.PredictOptions(cache=False, clinvar=False, protein_context=False, sequences=False)
    res = P.predict(text="7 140753336 A T", options=opts)
    assert set(res.table["nmd"]) == {"NMD uncertain (the change spans an exon-intron boundary)"}
