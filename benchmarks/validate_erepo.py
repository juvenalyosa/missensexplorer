"""Do the automated ACMG criteria agree with ClinGen's expert panels? A check on the ClinGen eRepo.

    python benchmarks/validate_erepo.py [--tsv erepo_classifications.tsv] [--limit N]

Data: the ClinGen Evidence Repository (eRepo; erepo.clinicalgenome.org, CC0), where ClinGen Variant
Curation Expert Panels (VCEPs) list, for every classified variant, the ACMG/AMP criteria they applied
("Met") and the ones they evaluated and found not met ("Not Met").

Plan (fixed before the results were seen; see the paper notes):
  * Non-retracted classifications of single-nucleotide missense variants (`c.N>N (p.Xxx#Yyy)`, not a
    stop and not synonymous), mapped to GRCh38 through their ClinVar Variation ID (variant_summary).
  * One classification per variant: the most recently approved.
  * MissensExplorer with its default ACMG settings: PP3/BP4 from the model score with its own
    ClinGen-method calibration; PS1/PM5 from live ClinVar, the variant itself left out; PM2_Supporting
    from a direct gnomAD check (confirmed absence and coverage); BA1 at gnomAD grpmax >= 5% with
    ClinGen's exception list.
  * Per criterion, among the classifications whose panel evaluated it (Met or Not Met):
    sensitivity (tool applies it when the panel did), specificity (tool does not apply it when the
    panel found it not met) and positive predictive value, with Wilson 95% intervals.
  * PP3/BP4: strength table (panel x tool) and opposite-direction conflicts.
  * Secondary: AUC of the model score, REVEL and AlphaMissense against the panel's final class
    (P/LP vs B/LB), split by overlap with the model's training data (variant in it / gene in it /
    gene never seen), with bootstrap intervals; and before 2023 vs 2023 or later (ClinGen's PP3/BP4
    calibration was published at the end of 2022).

Caveat: VCEPs apply gene-specific specifications (their own REVEL cut-offs, frequency thresholds for
PM2/BA1/BS1), so this measures agreement with expert practice, not correctness.
Writes benchmarks/results/erepo_validation.{md,json} and erepo_validation_variants.csv.gz; the
secondary analyses (per panel, clustered intervals, strength agreement, points) are in
erepo_analyses.py, which reads that table.

Changes after v2.4 (prompted by review): a variant's gene counts as "in training" when it is the gene
of any training variant's ClinVar name or any of its ClinVar gene symbols (an exact match on
'LDLR;MIR6886' missed LDLR); the model score is kept at full precision (it was the call's 2-decimal
score); proportions are rounded once, in the report.
"""

from __future__ import annotations

import argparse
import csv
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

from missensexplorer import PredictOptions, __version__, model, predict  # noqa: E402
from missensexplorer.cache import default_cache_dir  # noqa: E402
from missensexplorer.predict import _num  # noqa: E402

EREPO_URL = "https://erepo.clinicalgenome.org/evrepo/api/summary/classifications/download"
OUT = Path(__file__).resolve().parent / "results"
TRAINING = OUT / "model_training_set.csv.gz"
NAME = re.compile(r"^(NM_\d+\.\d+)\(([^)]+)\):c\.(\d+)([ACGT])>([ACGT]) \(p\.([A-Z][a-z]{2})(\d+)([A-Z][a-z]{2})\)$")
PATHOGENIC = {"Pathogenic", "Likely Pathogenic"}
BENIGN = {"Benign", "Likely Benign"}
CRITERIA = ("PP3", "BP4", "PM2", "PS1", "PM5", "BA1")
STRENGTH = {"Supporting": 1, "Moderate": 2, "Strong": 4, "Very Strong": 8, "Very": 8, "Stand Alone": 8}
DEFAULT_STRENGTH = {"PP3": 1, "BP4": 1, "PM2": 2, "PS1": 4, "PM5": 2, "BA1": 8}


def download(path: Path | None) -> Path:
    if path:
        return path
    f = default_cache_dir() / "benchmark" / "erepo_classifications.tsv"
    if not f.exists():
        f.parent.mkdir(parents=True, exist_ok=True)
        print("Downloading the ClinGen eRepo ...", flush=True)
        r = httpx.get(EREPO_URL, timeout=600, follow_redirects=True)
        r.raise_for_status()
        f.write_bytes(r.content)
    return f


def parse_codes(text: str) -> dict[str, int]:
    """'PP3_Moderate,PM2_Supporting,PS4' -> {'PP3': 2, 'PM2': 1, 'PS4': 4} (base code -> points)."""
    out: dict[str, int] = {}
    for c in (x.strip() for x in re.split(r"[,;]", text or "")):
        if not c:
            continue
        base, _, mod = c.partition("_")
        pts = STRENGTH.get(mod.strip(), DEFAULT_STRENGTH.get(base, 0)) if mod else DEFAULT_STRENGTH.get(base, 0)
        out[base] = max(out.get(base, 0), pts)
    return out


def load_erepo(path: Path) -> tuple[list[dict], dict]:
    csv.field_size_limit(10**9)
    with open(path, encoding="utf-8") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    steps = {"classifications": len(rows)}
    rows = [r for r in rows if r["Retracted"].strip().lower() != "true"]
    steps["not retracted"] = len(rows)
    keep = []
    for r in rows:
        m = NAME.match(r["Variation"].strip())
        if not m or m.group(8) in ("Ter", m.group(6)) or not r["ClinVar Variation Id"].strip().isdigit():
            continue
        keep.append(
            {
                "clinvar_id": r["ClinVar Variation Id"].strip(),
                "gene": r["HGNC Gene Symbol"].strip(),
                "name": r["Variation"].strip(),
                "assertion": r["Assertion"].strip(),
                "met": parse_codes(r["Applied Evidence Codes (Met)"]),
                "not_met": parse_codes(r["Applied Evidence Codes (Not Met)"]),
                "met_text": r["Applied Evidence Codes (Met)"].strip(),
                "panel": r["Expert Panel"].strip(),
                "approved": r["Approval Date"].strip(),
                "uuid": r["Uuid"].strip(),
            }
        )
    steps["single-nucleotide missense with a ClinVar Variation ID"] = len(keep)
    latest: dict[str, dict] = {}
    for k in keep:
        old = latest.get(k["clinvar_id"])
        if old is None or k["approved"] > old["approved"]:
            latest[k["clinvar_id"]] = k
    steps["unique variants (latest approval kept)"] = len(latest)
    return list(latest.values()), steps


def grch38_positions(ids: set[str]) -> dict[str, str]:
    out = {}
    with gzip.open(B.download(), "rt", encoding="utf-8", errors="replace") as fh:
        ix = {k: i for i, k in enumerate(fh.readline().lstrip("#").rstrip("\n").split("\t"))}
        for line in fh:
            f = line.rstrip("\n").split("\t")
            vid = f[ix["VariationID"]]
            if vid not in ids or f[ix["Assembly"]] != "GRCh38":
                continue
            chrom, pos = f[ix["Chromosome"]], f[ix["PositionVCF"]]
            ref, alt = f[ix["ReferenceAlleleVCF"]], f[ix["AlternateAlleleVCF"]]
            if chrom in B.CHROMS and pos.isdigit() and len(ref) == 1 and len(alt) == 1:
                out[vid] = f"{chrom}:{pos} {ref}>{alt}"
    return out


def tool_codes(r: dict) -> dict[str, int]:
    """The criteria the tool applied to one results row, as {base code: points}."""
    out: dict[str, int] = {}
    pts = _num(r.get("acmg_points"))
    if pts:
        out["PP3" if pts > 0 else "BP4"] = int(abs(pts))
    code = str(r.get("acmg_ps1_pm5") or "").split(" ")[0]
    if code.startswith(("PS1", "PM5")):
        base, _, mod = code.partition("_")
        out[base] = STRENGTH.get(mod, DEFAULT_STRENGTH[base]) if mod else DEFAULT_STRENGTH[base]
    if str(r.get("acmg_pm2") or "").startswith("PM2_Supporting"):
        out["PM2"] = 1
    if str(r.get("population_frequency") or "").startswith("common: BA1"):
        out["BA1"] = 8
    return out


def wilson(k: int, n: int) -> list[float | None]:
    if n == 0:
        return [None, None]
    z, p = 1.959963984540054, k / n
    d = 1 + z * z / n
    c = (p + z * z / (2 * n)) / d
    h = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d
    return [round(c - h, 6), round(c + h, 6)]


def agreement(df: pd.DataFrame, code: str) -> dict:
    ev = df[df[f"exp_{code}"].notna()]
    exp = ev[f"exp_{code}"] > 0
    tool = ev[f"tool_{code}"] > 0
    tp, fn = int((exp & tool).sum()), int((exp & ~tool).sum())
    fp, tn = int((~exp & tool).sum()), int((~exp & ~tool).sum())
    return {
        "evaluated": len(ev),
        "panel_met": tp + fn,
        "panel_not_met": fp + tn,
        "tp": tp,
        "fn": fn,
        "fp": fp,
        "tn": tn,
        "sensitivity": round(tp / (tp + fn), 6) if tp + fn else None,
        "sensitivity_ci": wilson(tp, tp + fn),
        "specificity": round(tn / (tn + fp), 6) if tn + fp else None,
        "specificity_ci": wilson(tn, tn + fp),
        "ppv": round(tp / (tp + fp), 6) if tp + fp else None,
        "ppv_ci": wilson(tp, tp + fp),
        "tool_applied_not_evaluated": int(((df[f"exp_{code}"].isna()) & (df[f"tool_{code}"] > 0)).sum()),
    }


def auc(score: np.ndarray, y: np.ndarray) -> float | None:
    ok = ~np.isnan(score)
    s, t = score[ok], y[ok]
    n1, n0 = int(t.sum()), int((t == 0).sum())
    if n1 == 0 or n0 == 0:
        return None
    ranks = pd.Series(s).rank().to_numpy()
    return float((ranks[t == 1].sum() - n1 * (n1 + 1) / 2) / (n1 * n0))


def auc_ci(score: np.ndarray, y: np.ndarray, rng: np.random.Generator, n_boot: int = 2000) -> dict:
    ok = ~np.isnan(score)
    s, t = score[ok], y[ok]
    est = auc(s, t)
    bs = []
    for _ in range(n_boot):
        i = rng.integers(0, len(t), len(t))
        a = auc(s[i], t[i])
        if a is not None:
            bs.append(a)
    lo, hi = np.percentile(bs, [2.5, 97.5]) if bs else (np.nan, np.nan)
    return {
        "n": len(t),
        "pathogenic": int(t.sum()),
        "benign": int((t == 0).sum()),
        "coverage": round(float(ok.mean()), 4) if len(ok) else None,
        "auc": None if est is None else round(est, 4),
        "ci": [round(float(lo), 4), round(float(hi), 4)],
    }


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--tsv", type=Path, help="eRepo tab-delimited export (default: download it)")
    ap.add_argument("--limit", type=int, help="first N variants only (for a quick test)")
    ap.add_argument("--workers", type=int, default=6)
    args = ap.parse_args(argv)
    t0 = time.time()
    path = download(args.tsv)
    sha = hashlib.sha256(path.read_bytes()).hexdigest()
    variants, steps = load_erepo(path)
    pos = grch38_positions({v["clinvar_id"] for v in variants})
    variants = [{**v, "input": pos[v["clinvar_id"]]} for v in variants if v["clinvar_id"] in pos]
    steps["mapped to GRCh38 (ClinVar variant_summary)"] = len(variants)
    variants.sort(key=lambda v: (v["gene"], v["clinvar_id"]))
    if args.limit:
        variants = variants[: args.limit]
    print(steps, flush=True)

    res = predict(
        text="\n".join(v["input"] for v in variants),
        options=PredictOptions(protein_context=False, sequences=False, workers=args.workers),
        progress=lambda s, d, n: print(f"  {s}: {d}/{n}      ", end="\r", flush=True),
    )
    print()
    failed = {k: v for k, v in res.sources.items() if str(v).startswith("failed")}
    if failed:
        raise SystemExit(f"annotation incomplete, not reporting it: {failed}. Run again (cached answers are kept).")

    by_input = defaultdict(list)
    for r in res.table.to_dict("records"):
        by_input[r["input"]].append(r)
    train = pd.read_csv(TRAINING, dtype={"clinvar_id": str})
    train_ids = set(train["clinvar_id"])
    train_genes = {g.strip() for s in train["gene"].astype(str) for g in s.split(";")}
    if "gene_key" in train:
        train_genes |= set(train["gene_key"].astype(str))
    recs = []
    for v in variants:
        cand = [r for r in by_input.get(v["input"], []) if "missense_variant" in str(r["consequence"])]
        same = [r for r in cand if r["gene"] == v["gene"]]
        r = (same or [None])[0]
        if r is None:
            continue
        tc = tool_codes(r)
        ms = model.score(r)
        rec = {
            "clinvar_id": v["clinvar_id"],
            "input": v["input"],
            "gene": v["gene"],
            "name": v["name"],
            "panel": v["panel"],
            "approved": v["approved"],
            "assertion": v["assertion"],
            "panel_codes": v["met_text"],
            "panel_met_points": json.dumps(v["met"], sort_keys=True),
            "panel_not_met": ",".join(sorted(v["not_met"])),
            "tool_codes": r.get("acmg_codes"),
            "tool_pp3_bp4": r.get("acmg_pp3_bp4"),
            "tool_ps1_pm5": r.get("acmg_ps1_pm5"),
            "tool_pm2": r.get("acmg_pm2"),
            "population_frequency": r.get("population_frequency"),
            "prediction": r.get("prediction"),
            "model_score": ms.p if ms is not None else None,
            "spliceai_max": _num(r.get("spliceai_max")),
            "revel": _num(r.get("revel")),
            "alphamissense": _num(r.get("alphamissense")),
            "grpmax": _num(r.get("gnomad_popmax_af")),
            "overlap": "variant in training"
            if v["clinvar_id"] in train_ids
            else ("gene in training" if v["gene"] in train_genes else "gene never seen"),
        }
        for c in CRITERIA:
            rec[f"exp_{c}"] = v["met"].get(c) if c in v["met"] else (0 if c in v["not_met"] else None)
            rec[f"tool_{c}"] = tc.get(c, 0)
        recs.append(rec)
    df = pd.DataFrame(recs)
    steps["annotated as missense on the panel's gene"] = len(df)
    print(steps, flush=True)

    rng = np.random.default_rng(0)
    out: dict = {
        "generated": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "missensexplorer": __version__,
        "run": {k: str(v) for k, v in (res.run or {}).items()} if isinstance(res.run, dict) else str(res.run),
        "erepo_sha256": sha,
        "steps": steps,
        "panels": int(df["panel"].nunique()),
        "genes": int(df["gene"].nunique()),
        "assertions": dict(Counter(df["assertion"])),
        "overlap": dict(Counter(df["overlap"])),
        "criteria": {},
        "criteria_by_period": {},
        "strength": {},
        "conflicts": {},
        "auc": {},
        "calls": {},
    }
    for c in CRITERIA:
        out["criteria"][c] = agreement(df, c)
    early = df["approved"] < "2023-01-01"
    for name, d in (("approved before 2023", df[early]), ("approved 2023 or later", df[~early])):
        out["criteria_by_period"][name] = {c: agreement(d, c) for c in ("PP3", "BP4", "PM2", "PS1", "PM5", "BA1")}
    for c in ("PP3", "BP4"):
        ev = df[df[f"exp_{c}"].notna()]
        out["strength"][c] = {
            f"panel {int(a)} / tool {int(b)}": int(n)
            for (a, b), n in ev.groupby([f"exp_{c}", f"tool_{c}"]).size().items()
        }
    out["conflicts"] = {
        "panel PP3, tool BP4": int(((df["exp_PP3"] > 0) & (df["tool_BP4"] > 0)).sum()),
        "panel BP4, tool PP3": int(((df["exp_BP4"] > 0) & (df["tool_PP3"] > 0)).sum()),
        "panel PP3 (any)": int((df["exp_PP3"] > 0).sum()),
        "panel BP4 (any)": int((df["exp_BP4"] > 0).sum()),
        "panel P/LP, tool BP4": int((df["assertion"].isin(PATHOGENIC) & (df["tool_BP4"] > 0)).sum()),
        "panel B/LB, tool PP3": int((df["assertion"].isin(BENIGN) & (df["tool_PP3"] > 0)).sum()),
        "panel P/LP": int(df["assertion"].isin(PATHOGENIC).sum()),
        "panel B/LB": int(df["assertion"].isin(BENIGN).sum()),
    }
    lab = df[df["assertion"].isin(PATHOGENIC | BENIGN)].copy()
    lab["y"] = lab["assertion"].isin(PATHOGENIC).astype(int)
    for name, d in (
        ("all", lab),
        ("variant not in training", lab[lab["overlap"] != "variant in training"]),
        ("gene never seen", lab[lab["overlap"] == "gene never seen"]),
        (
            "approved 2023 or later, variant not in training",
            lab[(lab["approved"] >= "2023-01-01") & (lab["overlap"] != "variant in training")],
        ),
    ):
        y = d["y"].to_numpy()
        out["auc"][name] = {
            s: auc_ci(d[s].astype(float).to_numpy(), y, rng) for s in ("model_score", "revel", "alphamissense")
        }
        out["calls"][name] = {
            cls: dict(Counter(d[d["assertion"].isin(grp)]["prediction"].fillna("none")))
            for cls, grp in (("P/LP", PATHOGENIC), ("B/LB", BENIGN))
        }
    vus = df[df["assertion"] == "Uncertain Significance"]
    out["calls"]["VUS (all)"] = dict(Counter(vus["prediction"].fillna("none")))
    out["seconds"] = round(time.time() - t0)

    OUT.mkdir(exist_ok=True)
    df.to_csv(OUT / "erepo_validation_variants.csv.gz", index=False, compression="gzip")
    (OUT / "erepo_validation.json").write_text(json.dumps(out, indent=1, default=str))
    (OUT / "erepo_validation.md").write_text(markdown(out))
    print(markdown(out))
    return 0


def _pct(x) -> str:
    return "—" if x is None else f"{100 * x:.1f}%"


def _ci(ci) -> str:
    return "" if ci[0] is None else f" ({100 * ci[0]:.1f}–{100 * ci[1]:.1f})"


def markdown(rep: dict) -> str:
    L = [
        "# Automated ACMG criteria against ClinGen expert panels (eRepo)",
        "",
        f"Generated by `benchmarks/validate_erepo.py` on {rep['generated']} with MissensExplorer "
        f"{rep['missensexplorer']}. eRepo export sha256 `{rep['erepo_sha256'][:16]}…`.",
        "",
        "Selection: " + "; ".join(f"{k}: {v:,}" for k, v in rep["steps"].items()) + ".",
        f"{rep['genes']} genes, {rep['panels']} expert panels. Panel classes: "
        + ", ".join(f"{k} {v:,}" for k, v in sorted(rep["assertions"].items()))
        + ". Overlap with the model's training data: "
        + ", ".join(f"{k} {v:,}" for k, v in sorted(rep["overlap"].items()))
        + ".",
        "",
        "Per criterion, among the classifications whose panel evaluated it (listed as met or not met). "
        "Positive agreement: the tool applies the criterion when the panel did. Negative agreement: the tool does not "
        "apply it when the panel found it not met. Wilson 95% intervals.",
        "",
        "| Criterion | Evaluated | Panel met | Positive agreement | Negative agreement | PPV "
        "| Tool applied, panel did not evaluate |",
        "|---|---|---|---|---|---|---|",
    ]
    for c, a in rep["criteria"].items():
        L.append(
            f"| {c} | {a['evaluated']:,} | {a['panel_met']:,} | {_pct(a['sensitivity'])}{_ci(a['sensitivity_ci'])} | "
            f"{_pct(a['specificity'])}{_ci(a['specificity_ci'])} | {_pct(a['ppv'])}{_ci(a['ppv_ci'])} | "
            f"{a['tool_applied_not_evaluated']:,} |"
        )
    for per, crit in rep["criteria_by_period"].items():
        L += [
            "",
            f"## {per}",
            "",
            "| Criterion | Evaluated | Panel met | Positive agreement | Negative agreement |",
            "|---|---|---|---|---|",
        ]
        for c, a in crit.items():
            L.append(
                f"| {c} | {a['evaluated']:,} | {a['panel_met']:,} | "
                f"{_pct(a['sensitivity'])}{_ci(a['sensitivity_ci'])} | "
                f"{_pct(a['specificity'])}{_ci(a['specificity_ci'])} |"
            )
    L += ["", "## PP3/BP4 strength (points: panel / tool)", ""]
    for c, tab in rep["strength"].items():
        L.append(f"- {c}: " + "; ".join(f"{k}: {v}" for k, v in tab.items()))
    L += ["", "## Opposite-direction evidence", ""]
    L += [f"- {k}: {v:,}" for k, v in rep["conflicts"].items()]
    L += [
        "",
        "## Ranking against the panel's class (P/LP vs B/LB)",
        "",
        "| Subset | Score | n (P/B) | Coverage | AUC (95% CI) |",
        "|---|---|---|---|---|",
    ]
    for sub, scores in rep["auc"].items():
        for s, a in scores.items():
            L.append(
                f"| {sub} | {s} | {a['n']} ({a['pathogenic']}/{a['benign']}) | {_pct(a['coverage'])} | "
                f"{a['auc']} ({a['ci'][0]}–{a['ci'][1]}) |"
            )
    L += ["", "## Model calls by panel class", ""]
    for sub, d in rep["calls"].items():
        L.append(f"- {sub}: {json.dumps(d)}")
    L += [
        "",
        "Caveats: VCEPs apply gene-specific specifications (own REVEL cut-offs; frequency thresholds for "
        "PM2/BA1/BS1 by disease), so the table measures agreement with expert practice, not correctness. "
        "PS1/PM5 use today's ClinVar, which can hold reference variants classified after the panel's "
        "decision. eRepo classifications are also in ClinVar (expert panel, 3 stars); variants in the model's "
        "training set are reported separately.",
        "",
    ]
    return "\n".join(L)


if __name__ == "__main__":
    sys.exit(main())
