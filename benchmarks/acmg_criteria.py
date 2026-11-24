"""How informative are the automated PS1, PM5 and PM2? A check on the held-out ClinVar benchmark.

    python benchmarks/acmg_criteria.py

The benchmark variants (clinvar_benchmark.py: 1,972 germline missense variants, >= 2 review stars,
from genes the model never saw) are annotated by the pipeline, including the direct gnomAD check
behind PM2. PS1/PM5 are computed by residue.py from ClinVar's variant_summary files instead of the
live E-utilities (the same records, much faster):
  * current ClinVar, the variant's own record left out, as in the app;
  * the January 2023 release, for the variants labelled since then (time split): their reference
    variants were classified before them.
For each criterion: how many pathogenic and benign variants get it, and its likelihood ratio (share of
pathogenic variants with it over share of benign ones) with a 95% bootstrap interval, against what its
strength needs in Bayes points (1 point 2.41, 2 points 5.79, 4 points 33.5; Pejaver et al. 2022).

Circularity: ClinVar submitters apply these same criteria (and frequency) when they classify, so the
ratios are optimistic; the time split limits this for PS1/PM5 but does not remove it.
Writes benchmarks/results/acmg_criteria.{md,json}.
"""

from __future__ import annotations

import gzip
import json
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(Path(__file__).resolve().parent))

import clinvar_benchmark as B  # noqa: E402
import compare_tools as C  # noqa: E402
import train_model as T  # noqa: E402

from missensexplorer import PredictOptions, predict, residue  # noqa: E402
from missensexplorer.cache import default_cache_dir  # noqa: E402
from missensexplorer.predict import _num  # noqa: E402

CODES = {"PS1": 4, "PS1_Moderate": 2, "PM5": 2, "PM5_Supporting": 1, "PM2_Supporting": 1}
NEEDED = {1: 2.41, 2: 5.79, 4: 33.53}


def pathogenic_index(path: Path) -> dict[tuple[str, int], dict[str, dict]]:
    """Index the reviewed pathogenic/likely pathogenic missense variants of a variant_summary file.

    Returns {(chrom, position): {VCV: record}} on GRCh38, the shape residue.references reads.
    """
    idx: dict[tuple[str, int], dict[str, dict]] = defaultdict(dict)
    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as fh:
        ix = {k: i for i, k in enumerate(fh.readline().lstrip("#").rstrip("\n").split("\t"))}
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if f[ix["Assembly"]] != "GRCh38" or f[ix["ClinicalSignificance"]] not in residue.PATHOGENIC:
                continue
            if f[ix["ReviewStatus"]] not in residue.REVIEWED or not residue.TITLE.match(f[ix["Name"]]):
                continue
            chrom, start, stop = f[ix["Chromosome"]], f[ix["Start"]], f[ix["Stop"]]
            if not (start.isdigit() and stop.isdigit()):
                continue
            rec = {
                "title": f[ix["Name"]],
                "classification": f[ix["ClinicalSignificance"]],
                "review_status": f[ix["ReviewStatus"]],
                "loc": {"GRCh38": [chrom, int(start), int(stop)]},
            }
            acc = f"VCV{int(f[ix['VariationID']]):09d}"
            for p in range(int(start), int(stop) + 1):
                idx[(chrom, p)][acc] = rec
    return idx


def residue_code(row: dict, idx) -> str | None:
    w = residue.window(row)
    if w is None:
        return None
    recs = {}
    for p in range(w[1], w[2] + 1):
        recs.update(idx.get((w[0], p), {}))
    text, _ = residue.evidence(row, residue.references(row, recs, "GRCh38"))
    return (text or "").split(" ")[0]


def lr(has: np.ndarray, y: np.ndarray, rng: np.random.Generator, n_boot: int = 2000) -> dict:
    def ratio(h, yy):
        a, b = h[yy == 1].mean(), h[yy == 0].mean()
        return np.inf if b == 0 else a / b

    idx = [rng.integers(0, len(y), len(y)) for _ in range(n_boot)]
    bs = np.array([ratio(has[i], y[i]) for i in idx])
    r = lambda v: "∞" if np.isinf(v) else round(float(v), 1)  # noqa: E731
    return {
        "pathogenic": int((has & (y == 1)).sum()),
        "benign": int((has & (y == 0)).sum()),
        "pathogenic_share": round(float(has[y == 1].mean()), 4),
        "benign_share": round(float(has[y == 0].mean()), 4),
        "lr": r(ratio(has, y)),
        "ci": [r(np.percentile(bs, 2.5, method="nearest")), r(np.percentile(bs, 97.5, method="nearest"))],
    }


def main() -> int:
    cands = B.load_candidates(B.download())
    picked = B.sample(cands, 1000, 8, 42)
    res = predict(text="\n".join(p["input"] for p in picked), options=PredictOptions(**T.OPTIONS), progress=T._progress)
    print()
    T.check_sources(res)
    old = C.labelled_in_2023()
    print("indexing ClinVar (current and January 2023) ...", flush=True)
    now = pathogenic_index(B.download())
    then = pathogenic_index(default_cache_dir() / "benchmark" / "variant_summary_2023-01.txt.gz")
    rows = []
    for p, r in T._missense_rows(res, picked):
        rows.append(
            {
                "id": p["id"],
                "label": p["label"],
                "new": p["id"] not in old,
                "rare": (_num(r.get("gnomad_popmax_af")) or 0) < T.RARE,
                "residue_now": residue_code(r, now),
                "residue_2023": residue_code(r, then),
                "pm2": str(r.get("acmg_pm2") or "").split(" ")[0],
            }
        )
    df = pd.DataFrame(rows)
    rng = np.random.default_rng(0)
    out = {
        "variants": len(df),
        "pathogenic": int(df["label"].sum()),
        "benign": int((df["label"] == 0).sum()),
        "new_since_2023": int(df["new"].sum()),
        "results": {},
    }
    sets = {
        "all variants, current ClinVar": (df, "residue_now"),
        "rare variants, current ClinVar": (df[df["rare"]], "residue_now"),
        "labelled since 2023, ClinVar of January 2023": (df[df["new"]], "residue_2023"),
    }
    for name, (d, col) in sets.items():
        y = d["label"].to_numpy()
        res_ = {"variants": len(d), "pathogenic": int(y.sum()), "benign": int((y == 0).sum())}
        for code in ("PS1", "PS1_Moderate", "PM5", "PM5_Supporting"):
            res_[code] = lr((d[col] == code).to_numpy(), y, rng)
        res_["any PS1/PM5"] = lr(d[col].isin(["PS1", "PS1_Moderate", "PM5", "PM5_Supporting"]).to_numpy(), y, rng)
        res_["PM2_Supporting"] = lr((d["pm2"] == "PM2_Supporting").to_numpy(), y, rng)
        out["results"][name] = res_
        print(
            name,
            {k: (v["pathogenic"], v["benign"], v["lr"]) for k, v in res_.items() if isinstance(v, dict)},
            flush=True,
        )
    T.OUT.mkdir(exist_ok=True)
    (T.OUT / "acmg_criteria.json").write_text(json.dumps(out, indent=1))
    (T.OUT / "acmg_criteria.md").write_text(markdown(out))
    return 0


def markdown(rep: dict) -> str:
    L = [
        "# Automated ACMG PS1, PM5 and PM2 on the held-out ClinVar benchmark",
        "",
        f"Generated by `benchmarks/acmg_criteria.py`: {rep['variants']:,} germline missense variants "
        f"({rep['pathogenic']:,} pathogenic, {rep['benign']:,} benign) from genes the model never saw; "
        f"{rep['new_since_2023']} labelled since 2023.",
        "",
        "Likelihood ratio = share of pathogenic variants with the criterion / share of benign ones. "
        "ClinVar submitters use these criteria themselves, so the ratios are optimistic; the time split "
        "(reference variants from the January 2023 release, tested on variants labelled later) limits this "
        "for PS1/PM5.",
        "",
    ]
    for name, r in rep["results"].items():
        L += [
            f"## {name}",
            "",
            f"{r['variants']:,} variants ({r['pathogenic']:,} pathogenic, {r['benign']:,} benign).",
            "",
            "| Criterion | Pathogenic with it | Benign with it | Likelihood ratio (95% CI) | Needed |",
            "|---|---|---|---|---|",
        ]
        for code in ("PS1", "PS1_Moderate", "PM5", "PM5_Supporting", "any PS1/PM5", "PM2_Supporting"):
            x = r[code]
            need = NEEDED.get(CODES.get(code), "—")
            L.append(
                f"| {code} | {x['pathogenic']} ({x['pathogenic_share']:.1%}) | {x['benign']} "
                f"({x['benign_share']:.1%}) | {x['lr']} ({x['ci'][0]}–{x['ci'][1]}) | {need} |"
            )
        L.append("")
    return "\n".join(L)


if __name__ == "__main__":
    sys.exit(main())
