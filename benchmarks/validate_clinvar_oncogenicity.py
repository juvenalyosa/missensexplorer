"""Do the somatic oncogenicity points agree with ClinVar's oncogenicity classifications?

    python benchmarks/validate_clinvar_oncogenicity.py [--limit N]

Data: ClinVar variant_summary (GRCh38), the oncogenicity classification that ClinVar aggregates from
submissions made under the ClinGen/CGC/VICC standard (Horak et al. 2022). Kept: classifications with
review criteria (at least one star), Oncogenic / Likely oncogenic / Uncertain significance / Likely
benign / Benign, variants of at most 50 bp with VCF coordinates.

The app never counts a variant's own ClinVar oncogenicity classification as evidence for itself, so the
comparison is not circular in that respect; other variants' classifications at the same residue do count
(OS1/OM4), as the standard intends, and how much the results rest on them is reported. The app's
default settings are used (gnomAD checked directly for variants without an Ensembl frequency).

Reported: the tool's class against ClinVar's (confusion table, rows = ClinVar, columns = tool), the
share of ClinVar oncogenic/likely oncogenic variants that the tool calls oncogenic/likely oncogenic,
opposite-direction calls, by consequence (every kind, including those with nothing recovered), by gene
role and its source, by ClinVar class (oncogenic vs likely oncogenic) and by submitter (ClinVar's
submission_summary); the classes with and without OP1/SBP1 in activating genes; the recovered variants
that need OS1/OM4 (other ClinVar classifications) or CIViC; a list of well-known hotspot variants; how
often each criterion fires; and, as an exploratory measure, how well the points separate
oncogenic/likely oncogenic from uncertain variants (AUC; an uncertain classification is not a negative
label). Writes benchmarks/results/clinvar_oncogenicity_validation.{md,json} and a per-variant table.
"""

from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import math
import re
import sys
import time
from collections import Counter, defaultdict
from pathlib import Path

import httpx
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(Path(__file__).resolve().parent))

import clinvar_benchmark as B  # noqa: E402

from missensexplorer import PredictOptions, __version__, predict  # noqa: E402
from missensexplorer.cache import default_cache_dir  # noqa: E402
from missensexplorer.oncogenicity import tier as onco_tier  # noqa: E402

OUT = Path(__file__).resolve().parent / "results"
ONC = {
    "Oncogenic": "O/LO",
    "Likely oncogenic": "O/LO",
    "Oncogenic/Likely oncogenic": "O/LO",
    "Uncertain significance": "VUS",
    "Likely benign": "B/LB",
    "Benign": "B/LB",
}
REVIEWED = {"criteria provided, single submitter", "criteria provided, multiple submitters, no conflicts"}
TOOL = {"Oncogenic": "O/LO", "Likely oncogenic": "O/LO", "VUS": "VUS", "Likely benign": "B/LB", "Benign": "B/LB"}
NULL = ("stop_gained", "frameshift_variant", "splice_donor_variant", "splice_acceptor_variant", "start_lost")
INFRAME = ("inframe_deletion", "inframe_insertion", "protein_altering_variant", "stop_lost")
SUBMISSIONS = "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/submission_summary.txt.gz"
LANDMARKS = [  # well-known oncogenic hotspots, shown one by one
    ("BRAF", "p.Val600Glu"),
    ("KRAS", "p.Gly12Cys"),
    ("KRAS", "p.Gly12Asp"),
    ("KRAS", "p.Gly12Val"),
    ("KRAS", "p.Gln61Arg"),
    ("NRAS", "p.Gln61Arg"),
    ("PIK3CA", "p.His1047Arg"),
    ("PIK3CA", "p.Glu545Lys"),
    ("IDH1", "p.Arg132His"),
    ("EGFR", "p.Leu858Arg"),
    ("CTNNB1", "p.Ser33Phe"),
    ("CTNNB1", "p.Ser37Phe"),
    ("CTNNB1", "p.Ser45Phe"),
    ("SF3B1", "p.Lys700Glu"),
    ("U2AF1", "p.Ser34Phe"),
    ("JAK2", "p.Val617Phe"),
    ("TP53", "p.Arg175His"),
    ("TP53", "p.Arg273His"),
    ("H3-3A", "p.Lys28Met"),
]


def load(path: Path) -> list[dict]:
    out, seen = [], set()
    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as fh:
        ix = {k: i for i, k in enumerate(fh.readline().lstrip("#").rstrip("\n").split("\t"))}
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if f[ix["Assembly"]] != "GRCh38":
                continue
            cls = ONC.get(f[ix["Oncogenicity"]])
            if cls is None or f[ix["ReviewStatusOncogenicity"]] not in REVIEWED:
                continue
            chrom, pos = f[ix["Chromosome"]], f[ix["PositionVCF"]]
            ref, alt = f[ix["ReferenceAlleleVCF"]], f[ix["AlternateAlleleVCF"]]
            vid = f[ix["VariationID"]]
            if chrom not in B.CHROMS or not pos.isdigit() or ref in ("", "na", "-") or alt in ("", "na", "-"):
                continue
            if max(len(ref), len(alt)) > 51 or vid in seen:
                continue
            seen.add(vid)
            out.append(
                {
                    "clinvar_id": vid,
                    "gene": f[ix["GeneSymbol"]].split(";")[0],
                    "name": f[ix["Name"]],
                    "clinvar_class": f[ix["Oncogenicity"]],
                    "clinvar_group": cls,
                    "evaluated": f[ix["OncogenicityLastEvaluated"]],
                    "input": f"{chrom}:{pos} {ref}>{alt}",
                }
            )
    return out


REF_ID = re.compile(r"(?:OS1|OM4)\S* \([^)]*?VCV0*(\d+)")  # ClinVar ids of the OS1/OM4 reference variants


STANDARD = re.compile(r"the standard's rule for OS1 would give (-?\d+) points")


def _standard_points(r: dict) -> float | None:
    """The points under the standard's OS1 rule (the tool keeps the stronger of OS1 and the same-change evidence)."""
    pts = r.get("oncogenicity_points")
    if not isinstance(pts, (int, float)):
        return None
    m = STANDARD.search(str(r.get("oncogenicity_evidence") or ""))
    return float(m.group(1)) if m else float(pts)


def sha256(path: Path) -> str | None:
    """SHA-256 of a file (None when it is missing)."""
    if not path.exists():
        return None
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for block in iter(lambda: fh.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def role_of(evidence: str) -> str:
    """The gene role shown in the evidence: 'loss-of-function (UniProt)', 'activating (IntOGen)', ...

    An IntOGen role left open because it rests on too few samples is 'undecided (IntOGen)', apart from the
    roles that IntOGen's cohorts leave ambiguous.
    """
    m = re.search(r"(activating|loss-of-function|ambiguous) \((UniProt|IntOGen)([^)]*)\)", evidence)
    if m:
        if m.group(1) == "ambiguous" and "too few to decide" in m.group(3):
            return f"undecided ({m.group(2)})"
        return f"{m.group(1)} ({m.group(2)})"
    return "not a known cancer gene" if "not a known cancer gene" in evidence else "unknown"


def kind(consequence: str) -> str:
    c = str(consequence)
    if any(k in c for k in NULL):
        return "null variant"
    if any(k in c for k in INFRAME):
        return "in-frame"
    if "missense_variant" in c:
        return "missense"
    return "other"


def submitters() -> dict[str, list[str]]:
    """VariationID -> submitters of its oncogenicity classifications (ClinVar submission_summary, cached)."""
    f = default_cache_dir() / "benchmark" / "submission_summary.txt.gz"
    if not f.exists():
        print("Downloading ClinVar's submission_summary (about 400 MB) ...", flush=True)
        with (
            httpx.stream("GET", SUBMISSIONS, timeout=1200, follow_redirects=True) as r,
            open(f.with_suffix(".part"), "wb") as fh,
        ):
            r.raise_for_status()
            for chunk in r.iter_bytes(1 << 20):
                fh.write(chunk)
        f.with_suffix(".part").rename(f)
    out: dict[str, list[str]] = defaultdict(list)
    with gzip.open(f, "rt", encoding="utf-8", errors="replace") as fh:
        header = None
        for line in fh:
            if line.startswith("##"):
                continue
            if line.startswith("#"):
                header = {k: i for i, k in enumerate(line.lstrip("#").rstrip("\n").split("\t"))}
                continue
            x = line.rstrip("\n").split("\t")
            onc = x[header["Oncogenicity"]]
            if onc not in ("-", "", "na"):
                out[x[header["VariationID"]]].append(x[header["Submitter"]])
    return dict(out)


def wilson(k: int, n: int) -> list[float | None]:
    if n == 0:
        return [None, None]
    z, p = 1.959963984540054, k / n
    d = 1 + z * z / n
    c = (p + z * z / (2 * n)) / d
    h = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d
    return [round(c - h, 6), round(c + h, 6)]


def auc(s: np.ndarray, y: np.ndarray) -> float | None:
    n1, n0 = int(y.sum()), int((y == 0).sum())
    if not n1 or not n0:
        return None
    r = pd.Series(s).rank().to_numpy()
    return float((r[y == 1].sum() - n1 * (n1 + 1) / 2) / (n1 * n0))


def tier(points: float) -> str:
    return "O/LO" if points >= 6 else "VUS" if points >= 0 else "B/LB"


def gene_ci(olo: pd.DataFrame, col: str, reps: int = 10000) -> list[float | None]:
    """95% interval of the share recovered from a bootstrap of whole genes (variants of a gene are not independent)."""
    if olo.empty:
        return [None, None]
    g = olo.assign(hit=(olo[col] == "O/LO").astype(float)).groupby("gene")["hit"].agg(["sum", "count"]).to_numpy()
    picks = np.random.default_rng(7).integers(0, len(g), (reps, len(g)))
    tot = g[picks].sum(1)
    share = tot[:, 0] / tot[:, 1]
    return [round(float(x), 6) for x in np.percentile(share, [2.5, 97.5])]


def summary(d: pd.DataFrame, col: str = "tool_group") -> dict:
    tab = pd.crosstab(d["clinvar_group"], d[col])
    olo = d[d["clinvar_group"] == "O/LO"]
    k = int((olo[col] == "O/LO").sum())
    b = d[d["clinvar_group"] == "B/LB"]
    return {
        "n": len(d),
        "by_clinvar": dict(Counter(d["clinvar_group"])),
        "confusion (rows ClinVar, columns tool)": {r: {c: int(n) for c, n in v.items()} for r, v in tab.iterrows()},
        "agreement_3class": round(float((d["clinvar_group"] == d[col]).mean()), 6) if len(d) else None,
        "olo_called_olo": [k, len(olo), round(k / len(olo), 6) if len(olo) else None, wilson(k, len(olo))],
        "olo_called_olo_gene_ci": gene_ci(olo, col),
        "olo_called_blb": int((olo[col] == "B/LB").sum()),
        "blb_called_olo": int((b[col] == "O/LO").sum()),
        "vus_called_olo": int(((d["clinvar_group"] == "VUS") & (d[col] == "O/LO")).sum()),
    }


def _codes(text: object) -> list[str]:
    return [x.strip().split(" ")[0] for x in str(text or "").split(",") if x.strip()]


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--limit", type=int)
    ap.add_argument("--workers", type=int, default=6)
    args = ap.parse_args(argv)
    t0 = time.time()
    path = B.download()
    variants = load(path)
    if args.limit:
        variants = variants[: args.limit]
    print(
        f"{len(variants)} ClinVar oncogenicity classifications with criteria",
        dict(Counter(v["clinvar_group"] for v in variants)),
        flush=True,
    )
    subs = submitters()
    res = predict(
        text="\n".join(v["input"] for v in variants),
        options=PredictOptions(oncogenicity=True, protein_context=False, sequences=False, workers=args.workers),
        progress=lambda s, d, n: print(f"  {s}: {d}/{n}      ", end="\r", flush=True),
    )
    print()
    failed = {k: v for k, v in res.sources.items() if str(v).startswith("failed")}
    if failed:
        raise SystemExit(f"annotation incomplete, not reporting it: {failed}. Run again (cached answers are kept).")
    by_input = defaultdict(list)
    for r in res.table.to_dict("records"):
        by_input[r["input"]].append(r)
    # the same variants without ClinVar's oncogenicity search (no OS1/OM4): the exact counterfactual
    print("again without ClinVar (OS1/OM4 counterfactual) ...", flush=True)
    off = predict(
        text="\n".join(v["input"] for v in variants),
        options=PredictOptions(
            oncogenicity=True, clinvar=False, protein_context=False, sequences=False, workers=args.workers
        ),
        progress=lambda s, d, n: print(f"  {s}: {d}/{n}      ", end="\r", flush=True),
    )
    print()
    off_points = {}
    for r in off.table.to_dict("records"):
        off_points.setdefault((r["input"], r.get("gene")), r.get("oncogenicity_points"))
    recs = []
    for v in variants:
        rows = by_input.get(v["input"], [])
        same = [r for r in rows if r.get("gene") == v["gene"]]
        r = (same or [None])[0]
        if r is None or not isinstance(r.get("oncogenicity"), str) or r["oncogenicity"].startswith("not assessed"):
            continue
        tclass = r["oncogenicity"]
        s = sorted(set(subs.get(v["clinvar_id"], [])))
        recs.append(
            {
                **{
                    k: v[k]
                    for k in ("clinvar_id", "gene", "name", "clinvar_class", "clinvar_group", "evaluated", "input")
                },
                "consequence": r.get("consequence"),
                "kind": kind(r.get("consequence")),
                "tool_class": tclass,
                "tool_group": TOOL.get(tclass, "other"),
                "points": r.get("oncogenicity_points"),
                "codes": r.get("oncogenicity_codes"),
                "own_clinvar_shown": r.get("clinvar_oncogenicity"),
                "reference_ids": ";".join(REF_ID.findall(str(r.get("oncogenicity_evidence") or ""))),
                "points_without_clinvar": off_points.get((v["input"], v["gene"])),
                "points_standard_os1": _standard_points(r),
                "evidence": r.get("oncogenicity_evidence"),
                "role": role_of(str(r.get("oncogenicity_evidence") or "")),
                "submitters": "; ".join(s),
            }
        )
    df = pd.DataFrame(recs)
    df["points"] = pd.to_numeric(df["points"], errors="coerce")
    # the classes without OP1/SBP1 in activating genes (the v2.4 design), for comparison
    act = df["role"].str.startswith("activating")
    adj = df["codes"].map(_codes).map(lambda c: ("OP1" in c) - ("SBP1" in c))
    df["points_no_op1_in_act"] = df["points"] - np.where(act, adj, 0)
    df["tool_group_no_op1_in_act"] = df["points_no_op1_in_act"].map(tier)
    # the variant's submitter with most oncogenicity classifications in this set (not the first by name)
    n_by = Counter(x for s in df["submitters"] for x in s.split("; ") if x)
    df["main_submitter"] = df["submitters"].map(
        lambda s: max((x for x in s.split("; ") if x), key=lambda x: (n_by[x], x), default="unknown")
    )
    out: dict = {
        "generated": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "missensexplorer": __version__,
        "clinvar_file": Path(path).name,  # the name only: the SHA-256 below identifies the file
        "selected": len(variants),
        "assessed": len(df),
        "genes": int(df["gene"].nunique()),
        "overall": summary(df),
        "by_kind": {k: summary(d) for k, d in df.groupby("kind")},
        "by_role": {k: summary(d) for k, d in df.groupby("role")},
        "by_clinvar_class": {k: summary(d) for k, d in df.groupby("clinvar_class")},
        "by_submitter": {k: summary(d) for k, d in df.groupby("main_submitter") if len(d) >= 20},
        "submitters_per_variant": dict(
            Counter(df["submitters"].str.count("; ") + (df["submitters"] != "").astype(int))
        ),
        "without_op1_in_activating_genes": {
            "overall": summary(df, "tool_group_no_op1_in_act"),
            "by_kind": {k: summary(d, "tool_group_no_op1_in_act") for k, d in df.groupby("kind")},
            "by_role": {k: summary(d, "tool_group_no_op1_in_act") for k, d in df.groupby("role")},
        },
        "criteria": {},
        "auc_olo_vs_vus": {},
    }
    rec = df[(df["clinvar_group"] == "O/LO") & (df["tool_group"] == "O/LO")]
    rc = rec["codes"].map(_codes)
    # exact: the points of the same variant assessed without ClinVar's oncogenicity search
    no_ref = pd.to_numeric(rec["points_without_clinvar"], errors="coerce")
    no_civic = rec["points"] - rc.map(lambda c: 4 * ("OS2" in c))
    out["recovered_dependence"] = {
        "recovered": len(rec),
        "with_os1_or_om4": int(rc.map(lambda c: any(x.startswith(("OS1", "OM4")) for x in c)).sum()),
        "still_recovered_without_os1_om4": int((no_ref >= 6).sum()),
        "with_os2": int(rc.map(lambda c: "OS2" in c).sum()),
        "still_recovered_without_os2": int((no_civic >= 6).sum()),
        "by_kind_needing_os1_om4": {k: int(((no_ref < 6) & (rec["kind"] == k)).sum()) for k in rec["kind"].unique()},
    }
    # OS1/OM4 references that are themselves in this evaluation set, classified by the same submitter
    ids = dict(zip(df["clinvar_id"].astype(str), df["submitters"], strict=True))
    used = rec[rc.map(lambda c: any(x.startswith(("OS1", "OM4")) for x in c))]
    in_set = used["reference_ids"].map(lambda r: [x for x in str(r).split(";") if x in ids])
    same_lab = [
        any(set(ids[x].split("; ")) & set(str(s).split("; ")) for x in refs)
        for refs, s in zip(in_set, used["submitters"], strict=True)
    ]
    needed = (no_ref < 6).reindex(used.index)
    out["recovered_dependence"]["os1_om4_reference_in_this_set"] = int(in_set.map(bool).sum())
    out["recovered_dependence"]["os1_om4_reference_in_this_set_same_submitter"] = int(sum(same_lab))
    out["recovered_dependence"]["decisive_with_same_submitter_reference"] = int(
        sum(bool(a) and bool(b) for a, b in zip(same_lab, needed, strict=True))
    )
    std = pd.to_numeric(df["points_standard_os1"], errors="coerce")
    std_group = std.map(lambda p: None if p != p else tier(p))
    differs = std != df["points"]
    out["standard_os1_rule"] = {
        "rows_with_other_points": int(differs.sum()),
        "class_changes": dict(
            Counter(
                f"{a} -> {b}" for a, b in zip(df.loc[differs, "tool_class"], std[differs].map(onco_tier), strict=True)
            )
        ),
        "oncogenic_calls": {
            "tool": int((df["tool_class"] == "Oncogenic").sum()),
            "standard_os1_rule": int((std.map(onco_tier) == "Oncogenic").sum()),
        },
        "olo_recovered": {
            "tool": int(((df["clinvar_group"] == "O/LO") & (df["tool_group"] == "O/LO")).sum()),
            "standard_os1_rule": int(((df["clinvar_group"] == "O/LO") & (std_group == "O/LO")).sum()),
        },
    }
    out["data_versions"] = {**res.run.get("data_versions", {}), "sources": res.sources}
    # the full record of both runs: software, Ensembl release, options, sources and data versions
    out["runs"] = {"with ClinVar": res.run, "without ClinVar": off.run}
    apc = df[(df["gene"] == "APC") & (df["kind"] == "null variant")]
    apc_codes = apc["codes"].map(_codes)
    out["apc_null_variants"] = {
        "n": len(apc),
        "ovs1": int(apc_codes.map(lambda c: "OVS1" in c).sum()),
        "ovs1_strong": int(apc_codes.map(lambda c: "OVS1_Strong" in c).sum()),
        "ovs1_moderate": int(apc_codes.map(lambda c: "OVS1_Moderate" in c).sum()),
        "recovered": int((apc["tool_group"] == "O/LO").sum()),
    }
    out["own_clinvar_classification_shown"] = int(df["own_clinvar_shown"].notna().sum())
    out["provenance"] = {
        "clinvar_file_sha256": sha256(Path(path)),
        "submission_summary_sha256": sha256(default_cache_dir() / "benchmark" / "submission_summary.txt.gz"),
    }
    out["landmarks"] = []
    for gene, prot in LANDMARKS:
        hit = df[(df["gene"] == gene) & df["name"].str.contains(re.escape(f"({prot})"))]
        for _, h in hit.head(2).iterrows():
            out["landmarks"].append(
                {
                    "variant": f"{gene} {prot}",
                    "clinvar": h["clinvar_class"],
                    "tool": h["tool_class"],
                    "points": h["points"],
                    "codes": h["codes"],
                }
            )
    codes = Counter(x for c in df["codes"].dropna() for x in _codes(c))
    out["criteria"] = dict(codes.most_common())
    sub = df[df["clinvar_group"].isin(["O/LO", "VUS"])].dropna(subset=["points"])
    for name, d in [("all", sub), *[(k, g) for k, g in sub.groupby("kind")]]:
        y = (d["clinvar_group"] == "O/LO").astype(int).to_numpy()
        a = auc(d["points"].to_numpy(float), y)
        out["auc_olo_vs_vus"][name] = {"n": len(d), "olo": int(y.sum()), "auc": None if a is None else round(a, 4)}
    out["seconds"] = round(time.time() - t0)
    OUT.mkdir(exist_ok=True)
    df.to_csv(OUT / "clinvar_oncogenicity_variants.csv.gz", index=False, compression="gzip")
    (OUT / "clinvar_oncogenicity_validation.json").write_text(json.dumps(out, indent=1, default=str))
    md = [
        "# Somatic oncogenicity points against ClinVar oncogenicity classifications",
        "",
        f"Generated by `benchmarks/validate_clinvar_oncogenicity.py` on {out['generated']} with MissensExplorer "
        f"{out['missensexplorer']}.",
        "",
        "```",
        json.dumps(out, indent=1, default=str),
        "```",
        "",
    ]
    (OUT / "clinvar_oncogenicity_validation.md").write_text("\n".join(md))
    print(
        json.dumps(
            {k: out[k] for k in ("selected", "assessed", "overall", "recovered_dependence")}, indent=1, default=str
        )
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
