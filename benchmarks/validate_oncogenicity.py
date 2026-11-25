"""Validation of the somatic oncogenicity points (oncogenicity.py) on laboratory data.

    python benchmarks/validate_oncogenicity.py [--bootstrap 1000]

Two genes with opposite roles, both measured in the laboratory:

- EGFR, an oncogene: Ikeuchi et al. (Oncogene 2022, CC BY 4.0) scored whether each EGFR variant seen
  in at least two tumours transforms NIH-3T3 cells (3-4 = transforming, 1 = not) and Ba/F3 cells
  (see validate_egfr.py).
- TP53, a tumour suppressor: the NCI TP53 Database (release 21) gives, for every TP53 missense SNV,
  the yeast transactivation class of Kato et al. (PNAS 2003; non-functional vs functional or
  supertrans) and the loss-of-function class from the human-cell screens of Giacomelli et al.
  (Nat Genet 2018).

Every variant goes through the full pipeline with ``oncogenicity=True``, as a user would run it, and
again without ClinVar (to see what the points do without ClinVar's oncogenicity classifications).
The points are compared with the model's score alone (AUC, ties count half), and the classes are
counted per assay class. TP53 is not in the model's training set; neither are the EGFR variants
(validate_egfr.py checks it). Only summaries are written; no per-variant table of the TP53 database.

CIViC evidence (OS2/SBS2 and the clinical summary) is used as the app uses it, except the items that
cite the assays themselves or yeast assays like Kato's (EXCLUDED_PMIDS): those would test the
points against their own source.
"""

from __future__ import annotations

import argparse
import json
import sys
from collections import Counter, defaultdict
from pathlib import Path

import httpx
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import validate_egfr as E

from missensexplorer import PredictOptions, model, predict
from missensexplorer import oncogenicity as onco
from missensexplorer.cache import default_cache_dir
from missensexplorer.oncogenicity import INTOGEN_RELEASE

TP53_DB = "https://tp53.cancer.gov/static/data/MutationView_r21.csv"
OUT = Path(__file__).resolve().parent / "results"
TIERS = ["Oncogenic", "Likely oncogenic", "VUS", "Likely benign", "Benign"]
EXCLUDED_PMIDS = {
    "35304574": "Ikeuchi et al. 2022 (the EGFR assay)",
    "12826609": "Kato et al. 2003 (the TP53 yeast assay)",
    "30224644": "Giacomelli et al. 2018 (the TP53 human-cell assay)",
    "10519380": "Marutani et al. 1999 (TP53 yeast assay)",
    "11896595": "Monti et al. 2002 (TP53 yeast assay)",
}
WATCH = {
    "EGFR": ["L858R", "G719A", "G719C", "G719S", "L861Q", "S768I", "E709K", "T790M", "L747P", "E746_A750del"],
    "TP53": ["R175H", "R248Q", "R248W", "R273H", "G245S", "Y220C", "R337H", "V217M", "P47S"],
}


def auc(y: np.ndarray, s: np.ndarray) -> float:
    """AUC with ties counted as half (Mann-Whitney), without sklearn."""
    order = np.argsort(s, kind="mergesort")
    ranks = np.empty(len(s))
    ss = s[order]
    i = 0
    while i < len(ss):
        j = i
        while j + 1 < len(ss) and ss[j + 1] == ss[i]:
            j += 1
        ranks[order[i : j + 1]] = (i + j) / 2 + 1
        i = j + 1
    npos = int(y.sum())
    nneg = len(y) - npos
    return float((ranks[y].sum() - npos * (npos + 1) / 2) / (npos * nneg))


def ci(y: np.ndarray, s: np.ndarray, rng: np.random.Generator, n: int) -> list[float]:
    vals = []
    for _ in range(n):
        i = rng.integers(0, len(y), len(y))
        if 0 < y[i].sum() < len(i):
            vals.append(auc(y[i], s[i]))
    return [round(float(np.percentile(vals, 2.5)), 3), round(float(np.percentile(vals, 97.5)), 3)]


_load = onco.load
dropped: set[str] = set()  # CIViC evidence IDs left out


def load_without_assay_papers(*args, **kwargs) -> onco.Data:
    """The app's data, without the CIViC items that cite the validation assays (see EXCLUDED_PMIDS)."""
    data = _load(*args, **kwargs)
    for gene, changes in data.civic.items():
        for change, items in changes.items():
            keep = [x for x in items if x.get("pmid") not in EXCLUDED_PMIDS]
            dropped.update(x["eid"] for x in items if x.get("pmid") in EXCLUDED_PMIDS)
            data.civic[gene][change] = keep
    return data


onco.load = load_without_assay_papers  # predict() reads it through the module


def run(inputs: list[str], clinvar: bool) -> pd.DataFrame:
    res = predict(
        text="\n".join(inputs),
        options=PredictOptions(
            oncogenicity=True, clinvar=clinvar, protein_context=False, sequences=False, gnomad_check=False, workers=6
        ),
        progress=lambda s, k, n: print(f"  {s}: {k}/{n}", end="\r", flush=True),
    )
    print()
    failed = {k: v for k, v in res.sources.items() if str(v).startswith("failed")}
    if failed:
        raise SystemExit(f"incomplete: {failed}. Run again (cached answers are kept).")
    return res.table


def summarise(rows: list[dict]) -> dict:
    """One record per assay variant: mean points and model score over its DNA changes, and its class."""
    if not rows:
        return {}
    pts = [r["oncogenicity_points"] for r in rows if isinstance(r.get("oncogenicity_points"), (int, float))]
    tiers = {r.get("oncogenicity") for r in rows}
    ps = [ms.p for r in rows if (ms := model.score(r)) is not None]
    return {
        "points": float(np.mean(pts)) if pts else None,
        "tier": tiers.pop() if len(tiers) == 1 else "mixed",
        "codes": " / ".join(sorted({str(r.get("oncogenicity_codes") or "") for r in rows})),
        "model": float(np.mean(ps)) if ps else None,
        "own_clinvar": any(isinstance(r.get("clinvar_oncogenicity"), str) for r in rows),
        "clinical_level": next(
            (r["clinical_evidence_level"] for r in rows if isinstance(r.get("clinical_evidence_level"), str)), None
        ),
    }


def egfr(client: httpx.Client) -> pd.DataFrame:
    path = default_cache_dir() / "validation" / "egfr_ikeuchi2022_table4.xlsx"
    path.parent.mkdir(parents=True, exist_ok=True)
    if not path.exists():
        r = client.get(E.TABLE)
        r.raise_for_status()
        path.write_bytes(r.content)
    items = []
    for x in E.read_table(path).to_dict("records"):
        kind, text = E.to_input(x["variant"])
        if text and kind in ("missense", "deletion", "insertion", "deletion-insertion"):
            items.append({"variant": x["variant"], "kind": kind, "input": text, "ffa": x["ffa"], "baf3": x["baf3"]})
    out = []
    for clinvar in (True, False):
        print(f"EGFR, ClinVar {'on' if clinvar else 'off'}: {len(items)} variants")
        by = defaultdict(list)
        for r in run(sorted({i["input"] for i in items}), clinvar).to_dict("records"):
            by[r["input"]].append(r)
        for i in items:
            s = summarise(E.pick_rows(i, by.get(i["input"], [])))
            out.append({**i, **s, "clinvar_on": clinvar})
    df = pd.DataFrame(out)
    df["y_3t3"] = df.ffa.map(lambda f: True if f >= 3 else False if f == 1 else None)
    df["y_baf3"] = df.baf3.astype(bool)
    df["assay_class"] = df.ffa.astype(str)
    return df


def tp53(client: httpx.Client) -> pd.DataFrame:
    path = default_cache_dir() / "validation" / "tp53_MutationView_r21.csv"
    if not path.exists():
        r = client.get(TP53_DB)
        r.raise_for_status()
        path.write_bytes(r.content)
    mv = pd.read_csv(path, low_memory=False)
    snv = mv[(mv.Effect == "missense") & mv.WT_nucleotide.str.len().eq(1) & mv.Mutant_nucleotide.str.len().eq(1)]
    snv = snv[snv.TransactivationClass.notna() | snv.DNE_LOFclass.isin(["DNE_LOF", "notDNE_LOF", "notDNE_notLOF"])]
    snv = snv.drop_duplicates("hgvs_hg38")
    out = []
    for clinvar in (True, False):
        print(f"TP53, ClinVar {'on' if clinvar else 'off'}: {len(snv)} SNVs")
        by = defaultdict(list)
        for r in run(list(snv.hgvs_hg38), clinvar).to_dict("records"):
            by[r["input"]].append(r)
        for prot, g in snv.groupby("ProtDescription"):
            v = prot.removeprefix("p.")
            rows = [
                r
                for h in g.hgvs_hg38
                for r in by.get(h, [])
                if r.get("gene") == "TP53"
                and r.get("protein_change") == v
                and "missense_variant" in str(r.get("consequence"))
            ]
            first = g.iloc[0]
            out.append(
                {
                    "variant": v,
                    "kind": "missense",
                    "kato": first.TransactivationClass if isinstance(first.TransactivationClass, str) else "",
                    "giacomelli": first.DNE_LOFclass if isinstance(first.DNE_LOFclass, str) else "",
                    **summarise(rows),
                    "clinvar_on": clinvar,
                }
            )
    df = pd.DataFrame(out)
    df["y_kato"] = df.kato.map({"non-functional": True, "functional": False, "supertrans": False})
    df["y_giacomelli"] = df.giacomelli.map({"DNE_LOF": True, "notDNE_LOF": True, "notDNE_notLOF": False})
    df["assay_class"] = df.kato.replace("", "no Kato class")
    return df


def evaluate(df: pd.DataFrame, labels: dict[str, str], classes: list[str], rng, n_boot: int) -> dict:
    on, off = df[df.clinvar_on], df[~df.clinvar_on].set_index("variant")
    rep: dict = {"variants": len(on), "assessed": int(on.points.notna().sum()), "auc": {}, "tiers": {}}
    for name, col in labels.items():
        for kind in ("missense", "all"):
            d = on[on[col].notna() & on.points.notna() & ((on.kind == "missense") if kind == "missense" else True)]
            if d.empty or (kind == "all" and (on.kind == "missense").all()):
                continue
            y = d[col].astype(bool).to_numpy()
            line = {"n": len(y), "positives": int(y.sum())}
            pts = d.points.to_numpy(float)
            line["points"] = round(auc(y, pts), 3)
            line["points_ci"] = ci(y, pts, rng, n_boot)
            off_pts = off.loc[d.variant, "points"].to_numpy(float)
            keep = ~np.isnan(off_pts)
            line["points_without_clinvar"] = round(auc(y[keep], off_pts[keep]), 3)
            m = d.model.notna().to_numpy()
            if m.sum() and 0 < y[m].sum() < m.sum():
                line["model"] = round(auc(y[m], d.model.to_numpy(float)[m]), 3)
                line["model_n"] = int(m.sum())
            rep["auc"][f"{name} ({kind})"] = line
    for c in classes:
        cnt = Counter(on.loc[on.assay_class == c, "tier"].fillna("not assessed"))
        rep["tiers"][c] = {t: int(cnt.get(t, 0)) for t in [*TIERS, "mixed", "not assessed"] if cnt.get(t)}
    codes = on.codes.fillna("")
    rep["civic"] = {
        "OS2": int(codes.str.contains(r"\bOS2\b").sum()),
        "SBS2": int(codes.str.contains(r"\bSBS2\b").sum()),
        "clinical_evidence": {
            lv: int((on.clinical_level == lv).sum()) for lv in "ABCDE" if (on.clinical_level == lv).any()
        },
    }
    return rep


def markdown(rep: dict, egfr_df: pd.DataFrame, tp53_df: pd.DataFrame) -> str:
    L = [
        "# Somatic oncogenicity points on laboratory data (EGFR and TP53)",
        "",
        "Produced by `benchmarks/validate_oncogenicity.py`. The points follow the ClinGen/CGC/VICC standard "
        "(Horak et al., Genet Med 2022) with the criteria computed automatically (see "
        "`src/missensexplorer/oncogenicity.py`): hotspots from cancerhotspots.org, gene roles from IntOGen "
        f"({INTOGEN_RELEASE}), ClinVar oncogenicity classifications of other variants at the residue, "
        "functional evidence curated in CIViC (OS2/SBS2), gnomAD, and the model's call (only in genes IntOGen "
        "does not call activating).",
        "",
        "- **EGFR** (oncogene): Ikeuchi et al., Oncogene 2022 (CC BY 4.0): NIH-3T3 transformation "
        "(3-4 vs 1) and Ba/F3 transformation.",
        "- **TP53** (tumour suppressor): NCI TP53 Database r21: Kato et al. 2003 yeast transactivation "
        "(non-functional vs functional/supertrans) and Giacomelli et al. 2018 loss of function.",
        "- **AUC**: chance that a variant positive in the assay has more points than a negative one (ties "
        "count half), with a 95% bootstrap interval. *Without ClinVar*: the same run with ClinVar off.",
        "- The points are an integer scale: many variants share a value, which caps the AUC.",
        "",
        "## Ranking",
        "",
        "| Gene | Assay | Variants | Positives | Model alone | Points | Points without ClinVar |",
        "|---|---|---|---|---|---|---|",
    ]
    for gene in ("EGFR", "TP53"):
        for k, v in rep[gene]["auc"].items():
            L.append(
                f"| {gene} | {k} | {v['n']} | {v['positives']} | {v.get('model', '–')} | "
                f"**{v['points']}** ({v['points_ci'][0]}–{v['points_ci'][1]}) | {v['points_without_clinvar']} |"
            )
    L += [
        "",
        "## CIViC",
        "",
        "Items that cite the validation assays or yeast assays like Kato's are left out here (not in the app): "
        + "; ".join(f"PMID {k}, {v}" for k, v in EXCLUDED_PMIDS.items())
        + f". Items left out: {rep['civic_dropped']}.",
        "",
        "| Gene | OS2 (functional evidence for) | SBS2 (against) | Clinical evidence, by best level |",
        "|---|---|---|---|",
    ]
    for gene in ("EGFR", "TP53"):
        c = rep[gene]["civic"]
        lv = ", ".join(f"{k} {n}" for k, n in c["clinical_evidence"].items()) or "none"
        L.append(f"| {gene} | {c['OS2']} | {c['SBS2']} | {lv} |")
    L += ["", "## Classes", ""]
    for gene, title in (("EGFR", "3T3 score"), ("TP53", "Kato class")):
        cols = [*TIERS, "mixed", "not assessed"]  # mixed: its DNA changes get different classes
        L += [f"**{gene}**, by {title}:", "", f"| {title} | " + " | ".join(cols) + " |", "|---" * (len(cols) + 1) + "|"]
        for c, cnt in rep[gene]["tiers"].items():
            L.append(f"| {c} | " + " | ".join(str(cnt.get(t, 0)) for t in cols) + " |")
        L.append("")
    L += ["## Examples", "", "| Gene | Variant | Assay | Class | Points | Criteria |", "|---|---|---|---|---|---|"]
    for gene, df, col in (("EGFR", egfr_df, "assay_class"), ("TP53", tp53_df, "kato")):
        d = df[df.clinvar_on].set_index("variant")
        for v in WATCH[gene]:
            if v in d.index:
                r = d.loc[v]
                assay = f"3T3 {r[col]}" if gene == "EGFR" else r[col]
                pts = "–" if pd.isna(r.points) else f"{r.points:g}"
                L.append(f"| {gene} | {v} | {assay} | {r.tier} | {pts} | {r.codes} |")
    L += [
        "",
        "Variants whose own allele has a ClinVar oncogenicity classification keep it in the ClinVar column; "
        "it never counts as evidence for itself. Insertions written as protein changes (most EGFR exon 20 "
        "insertions) cannot be traced to one DNA change and are not assessed.",
    ]
    return "\n".join(L) + "\n"


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--bootstrap", type=int, default=1000)
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--dump", type=Path, help="also write the per-variant tables here (for checking; not shared)")
    args = ap.parse_args(argv)
    rng = np.random.default_rng(args.seed)
    with httpx.Client(timeout=120, follow_redirects=True, headers={"User-Agent": "missensexplorer-benchmark"}) as c:
        e, t = egfr(c), tp53(c)
    rep = {
        "EGFR": evaluate(e, {"3T3": "y_3t3", "Ba/F3": "y_baf3"}, ["4", "3", "2", "1"], rng, args.bootstrap),
        "TP53": evaluate(
            t,
            {"Kato": "y_kato", "Giacomelli": "y_giacomelli"},
            ["non-functional", "partially functional", "functional", "supertrans"],
            rng,
            args.bootstrap,
        ),
        "sources": {"egfr": E.PAPER, "tp53": TP53_DB, "intogen": INTOGEN_RELEASE, "civic": onco.CIVIC_URL},
        "civic_dropped": len(dropped),
    }
    if args.dump:
        args.dump.mkdir(parents=True, exist_ok=True)
        e.to_csv(args.dump / "egfr.csv", index=False)
        t.to_csv(args.dump / "tp53.csv", index=False)
    OUT.mkdir(exist_ok=True)
    (OUT / "oncogenicity_validation.json").write_text(json.dumps(rep, indent=1))
    (OUT / "oncogenicity_validation.md").write_text(markdown(rep, e, t))
    print(json.dumps({g: rep[g]["auc"] for g in ("EGFR", "TP53")}, indent=1))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
