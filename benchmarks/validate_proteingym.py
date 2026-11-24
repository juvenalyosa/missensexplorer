"""Independent validation: MissensExplorer and other predictors against deep mutational scanning data.

    python benchmarks/validate_proteingym.py [--per-assay 200]

ClinVar labels are used to train and benchmark the model; this check uses none. ProteinGym
(Notin et al., NeurIPS 2023; https://proteingym.org) collects deep mutational scanning (DMS) assays:
laboratory measurements of how each amino-acid substitution changes a protein's function.

1. Downloads the ProteinGym substitution assays (v1.3, ~45 MB) and keeps the human ones.
2. Maps each assay onto the Ensembl canonical transcript of its gene (the assay's sequence must occur
   exactly in the canonical protein) and turns each substitution into a single-nucleotide variant
   using the coding sequence; substitutions that need two or three nucleotide changes are dropped,
   as the predictors only score single-nucleotide variants.
3. Annotates a random sample of each assay with the pipeline and with dbNSFP (as compare_tools.py).
4. Reports, per tool, the Spearman correlation between its score and the measured loss of function,
   averaged over the assays of each protein and then over proteins (ProteinGym's convention), with a
   95% bootstrap interval over proteins and the paired difference from MissensExplorer; and the AUC
   against ProteinGym's binary fit/unfit label.

DMS assays measure one molecular function in one experimental system, not disease, so correlations
are moderate for every tool; what matters is the ranking between tools on the same variants.
"""

from __future__ import annotations

import argparse
import io
import itertools
import json
import random
import sys
import time
import zipfile
from collections import defaultdict
from pathlib import Path

import httpx
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
sys.path.insert(0, str(Path(__file__).resolve().parent))

import compare_tools as C

from missensexplorer import PredictOptions, model, predict
from missensexplorer.cache import default_cache_dir
from missensexplorer.model_params import MODEL
from missensexplorer.predict import _num

PG = "https://marks.hms.harvard.edu/proteingym/ProteinGym_v1.3/DMS_ProteinGym_substitutions.zip"
PG_REF = "https://raw.githubusercontent.com/OATML-Markslab/ProteinGym/main/reference_files/DMS_substitutions.csv"
ENSEMBL = "https://rest.ensembl.org"
OUT = Path(__file__).resolve().parent / "results"
CODE = {
    "".join(c): aa
    for c, aa in zip(
        itertools.product("TCAG", repeat=3),
        "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
        strict=True,
    )
}
COMP = str.maketrans("ACGT", "TGCA")


def _get(client: httpx.Client, url: str, **kw):
    err = None
    for i in range(6):
        try:
            r = client.get(url, timeout=180, **kw)
            if r.status_code == 200:
                return r
            err = r.status_code
        except httpx.HTTPError as e:
            err = e
        time.sleep(2 * (i + 1))
    raise SystemExit(f"{url}: {err}")


def download(d: Path, client: httpx.Client) -> pd.DataFrame:
    d.mkdir(parents=True, exist_ok=True)
    if not (d / "DMS_substitutions.csv").exists():
        (d / "DMS_substitutions.csv").write_bytes(_get(client, PG_REF).content)
    if not (d / "proteingym").exists():
        print("Downloading ProteinGym substitutions ...", flush=True)
        zipfile.ZipFile(io.BytesIO(_get(client, PG, follow_redirects=True).content)).extractall(d / "proteingym")
    ref = pd.read_csv(d / "DMS_substitutions.csv")
    return ref[ref.taxon == "Human"]


def uniprot_genes(ref: pd.DataFrame, d: Path, client: httpx.Client) -> dict:
    path = d / "uniprot_genes.json"
    out = json.loads(path.read_text()) if path.exists() else {}
    for name in sorted(set(ref.UniProt_ID) - set(out)):
        res = (
            _get(
                client,
                "https://rest.uniprot.org/uniprotkb/search",
                params={"query": f"id:{name} AND reviewed:true", "fields": "accession,gene_primary", "format": "json"},
            )
            .json()
            .get("results", [])
        )
        out[name] = (res[0].get("genes") or [{}])[0].get("geneName", {}).get("value") if res else None
        path.write_text(json.dumps(out))
    return out


def transcript(gene: str, d: Path, client: httpx.Client) -> dict:
    """Canonical transcript of a gene: protein, CDS and the genomic position of every CDS base."""
    path = d / "transcripts" / f"{gene}.json"
    if path.exists():
        return json.loads(path.read_text())
    h = {"Content-Type": "application/json"}
    g = _get(client, f"{ENSEMBL}/lookup/symbol/homo_sapiens/{gene}", headers=h).json()
    tx = g["canonical_transcript"].split(".")[0]
    prot = _get(client, f"{ENSEMBL}/sequence/id/{tx}?type=protein", headers=h).json()["seq"]
    cds = _get(client, f"{ENSEMBL}/sequence/id/{tx}?type=cds", headers=h).json()["seq"]
    gpos: list[int | None] = []
    for m in _get(client, f"{ENSEMBL}/map/cds/{tx}/1..{len(cds)}", headers=h).json()["mappings"]:
        n = m["end"] - m["start"] + 1
        if m.get("gap"):
            gpos += [None] * n
        else:
            gpos += list(range(m["start"], m["end"] + 1) if m["strand"] == 1 else range(m["end"], m["start"] - 1, -1))
    t = {"tx": tx, "chrom": g["seq_region_name"], "strand": g["strand"], "prot": prot, "cds": cds, "gpos": gpos}
    path.parent.mkdir(exist_ok=True)
    path.write_text(json.dumps(t))
    return t


def snv_for(t: dict, p: int, wt: str, mt: str) -> str | None:
    """First single-nucleotide change of codon p that turns wt into mt, as 'chrom:pos REF>ALT'."""
    c0 = 3 * (p - 1)
    codon = t["cds"][c0 : c0 + 3]
    if CODE.get(codon) != wt:
        return None
    for i in range(3):
        for b in "ACGT":
            if b != codon[i] and CODE.get(codon[:i] + b + codon[i + 1 :]) == mt and t["gpos"][c0 + i]:
                rb, ab = (codon[i], b) if t["strand"] == 1 else (codon[i].translate(COMP), b.translate(COMP))
                return f"{t['chrom']}:{t['gpos'][c0 + i]} {rb}>{ab}"
    return None


def collect(ref: pd.DataFrame, genes: dict, d: Path, client: httpx.Client, per_assay: int, seed: int):
    rng = random.Random(seed)
    items, skipped = [], []
    for _, a in ref.iterrows():
        gene = genes.get(a.UniProt_ID)
        if not gene:
            skipped.append([a.DMS_id, "no gene"])
            continue
        t = transcript(gene, d, client)
        off = t["prot"].find(a.target_seq)
        if off < 0:
            skipped.append([a.DMS_id, "assay sequence not in the canonical protein"])
            continue
        rows = pd.read_csv(d / "proteingym" / "DMS_ProteinGym_substitutions" / a.DMS_filename)
        rows = rows[~rows.mutant.str.contains(":")].to_dict("records")
        rng.shuffle(rows)
        n = 0
        for m in rows:
            wt, pos, mt = m["mutant"][0], int(m["mutant"][1:-1]), m["mutant"][-1]
            if mt == "*" or a.target_seq[pos - 1] != wt:
                continue
            snv = snv_for(t, pos + off, wt, mt)
            if snv is None:
                continue
            items.append(
                {
                    "assay": a.DMS_id,
                    "protein": a.UniProt_ID,
                    "gene": gene,
                    "change": f"{wt}{pos + off}{mt}",
                    "snv": snv,
                    "dms": m["DMS_score"],
                    "bin": m["DMS_score_bin"],
                }
            )
            n += 1
            if n >= per_assay:
                break
    return items, skipped


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--per-assay", type=int, default=200)
    ap.add_argument("--seed", type=int, default=1)
    args = ap.parse_args(argv)
    d = default_cache_dir() / "validation"
    client = httpx.Client(headers={"User-Agent": "missensexplorer-benchmark"})
    ref = download(d, client)
    genes = uniprot_genes(ref, d, client)
    items, skipped = collect(ref, genes, d, client, args.per_assay, args.seed)
    inputs = sorted({i["snv"] for i in items})
    print(
        f"{len(items)} substitutions from {len({i['assay'] for i in items})} human assays; {len(inputs)} variants",
        flush=True,
    )

    res = predict(
        text="\n".join(inputs),
        options=PredictOptions(clinvar=False, protein_context=False, sequences=False, workers=6),
        progress=lambda s, k, n: print(f"  {s}: {k}/{n}", end="\r", flush=True),
    )
    print()
    failed = {k: v for k, v in res.sources.items() if v.startswith("failed")}
    if failed:
        raise SystemExit(f"annotation incomplete: {failed}. Run again (cached answers are kept).")
    by = defaultdict(list)
    for r in res.table.to_dict("records"):
        by[r["input"]].append(r)
    gene_change = {i["snv"]: (i["gene"], i["change"]) for i in items}
    ann = {}
    for inp in inputs:
        g, ch = gene_change[inp]
        cand = [
            r
            for r in by.get(inp, [])
            if "missense_variant" in str(r["consequence"]) and r["gene"] == g and r.get("protein_change") == ch
        ]
        if cand:
            ann[inp] = cand[0]
    ext = C.fetch_dbnsfp([{"id": k, "gene": gene_change[k][0], "input": k} for k in ann])

    names = ["MissensExplorer v" + MODEL["version"], *list(C.TOOLS)]
    per_assay = defaultdict(lambda: defaultdict(list))
    for it in items:
        r = ann.get(it["snv"])
        if r is None:
            continue
        ms = model.score(r)
        vals = [ms.p if ms else None]
        for name, (src, direction, _, _) in C.TOOLS.items():
            v = _num(r.get(src[4:])) if src.startswith("row:") else ext[it["snv"]].get(name)
            vals.append(None if v is None else direction * v)
        for n, v in zip(names, vals, strict=True):
            per_assay[it["assay"]][n].append((v, it["dms"], it["bin"], vals[0]))
    protein_of = {it["assay"]: it["protein"] for it in items}
    report = evaluate(names, per_assay, protein_of, np.random.default_rng(args.seed))
    report.update(
        {
            "proteingym": PG,
            "per_assay": args.per_assay,
            "assays": len(per_assay),
            "proteins": len({protein_of[a] for a in per_assay}),
            "variants": sum(len(v[names[0]]) for v in per_assay.values()),
            "skipped": skipped,
            "run": res.run,
        }
    )
    OUT.mkdir(exist_ok=True)
    (OUT / "proteingym_validation.json").write_text(json.dumps(report, indent=1, default=str))
    (OUT / "proteingym_validation.md").write_text(markdown(report))
    print((OUT / "proteingym_validation.md").read_text())
    return 0


def evaluate(names, per_assay, protein_of, rng) -> dict:
    from scipy.stats import spearmanr
    from sklearn.metrics import roc_auc_score

    per_protein = defaultdict(lambda: defaultdict(lambda: defaultdict(list)))  # tool -> protein -> metric -> values
    for assay, by_tool in per_assay.items():
        for n, vals in by_tool.items():
            ok = [(s, dms, b) for s, dms, b, ours in vals if s is not None and ours is not None]
            if len(ok) < 20:
                continue
            s, dms, b = (np.array(x, dtype=float) for x in zip(*ok, strict=True))
            rho = spearmanr(s, -dms).statistic  # higher score = more damaging; lower DMS score = less function
            if rho == rho:
                per_protein[n][protein_of[assay]]["rho"].append(rho)
            if (b == 0).sum() >= 5 and (b == 1).sum() >= 5:
                per_protein[n][protein_of[assay]]["auc"].append(roc_auc_score(1 - b, s))
    ours = {p: np.mean(m["rho"]) for p, m in per_protein[names[0]].items() if m["rho"]}
    tools = {}
    for n in names:
        prot = {p: np.mean(m["rho"]) for p, m in per_protein[n].items() if m["rho"]}
        common = sorted(set(prot) & set(ours))
        if len(common) < 10:
            continue
        a, o = np.array([prot[p] for p in common]), np.array([ours[p] for p in common])
        idx = [rng.integers(0, len(a), len(a)) for _ in range(2000)]
        aucs = [np.mean(per_protein[n][p]["auc"]) for p in common if per_protein[n][p]["auc"]]
        tools[n] = {
            "labels": "ClinVar (other genes)" if n == names[0] else C.TOOLS[n][2],
            "proteins": len(common),
            "spearman": round(float(a.mean()), 4),
            "ci": [round(float(np.percentile([a[i].mean() for i in idx], q)), 4) for q in (2.5, 97.5)],
            "difference": round(float((a - o).mean()), 4),
            "difference_ci": [
                round(float(np.percentile([(a[i] - o[i]).mean() for i in idx], q)), 4) for q in (2.5, 97.5)
            ],
            "auc": round(float(np.mean(aucs)), 4) if aucs else None,
        }
    return {"tools": dict(sorted(tools.items(), key=lambda kv: -kv[1]["spearman"]))}


def markdown(rep: dict) -> str:
    L = [
        "# Validation on deep mutational scanning data (ProteinGym)",
        "",
        f"{rep['variants']} single-nucleotide missense variants from {rep['assays']} human assays on "
        f"{rep['proteins']} proteins of [ProteinGym]({rep['proteingym']}) (up to {rep['per_assay']} per assay). "
        "No ClinVar label is involved: each variant's effect was measured in the laboratory. Produced by "
        "`benchmarks/validate_proteingym.py`.",
        "",
        "- **Spearman**: correlation between the tool's score and the measured loss of function, averaged over "
        "the assays of each protein and then over proteins. 95% bootstrap interval over proteins.",
        "- **Δ vs ours**: the tool's mean minus MissensExplorer's on the same proteins and variants; an interval "
        "that excludes 0 is a real difference.",
        "- **AUC**: against ProteinGym's binary label (fit/unfit), same averaging.",
        "- DMS assays measure one molecular function in one experimental system, not disease, so every tool's "
        "correlation is moderate. VARITY used some DMS data in training.",
        "",
        "| Tool | Clinical training labels | Proteins | Spearman | Δ vs ours | AUC |",
        "|---|---|---|---|---|---|",
    ]
    for n, t in rep["tools"].items():
        bold = "**" if n.startswith("MissensExplorer") else ""
        delta = "—" if bold else f"{t['difference']:+.3f} [{t['difference_ci'][0]:+.3f}, {t['difference_ci'][1]:+.3f}]"
        auc = "—" if t["auc"] is None else f"{t['auc']:.3f}"
        L.append(
            f"| {bold}{n}{bold} | {t['labels']} | {t['proteins']} | {t['spearman']:.3f} "
            f"[{t['ci'][0]:.3f}–{t['ci'][1]:.3f}] | {delta} | {auc} |"
        )
    L += [
        "",
        f"Assays left out ({len(rep['skipped'])}): their sequence does not occur exactly in the gene's "
        "canonical protein (isolated domains with tags, engineered constructs) — "
        + ", ".join(s[0] for s in rep["skipped"])
        + ".",
        "",
    ]
    return "\n".join(L)


if __name__ == "__main__":
    sys.exit(main())
