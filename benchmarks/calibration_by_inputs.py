"""PP3/BP4 band likelihood ratios by which predictor scores are missing (offline, cached data).

    python benchmarks/calibration_by_inputs.py

The model fills a missing score with its training mean (a zero term), so a band could be calibrated
overall but not for variants with every score, or not for those missing one. This reports the band
likelihood ratios of the shipped model on the ClinVar 2020 test set of Pejaver et al. (training variants
removed, as in calibrate_acmg.py), separately for variants with all four scores, with popEVE missing
(AlphaMissense and REVEL present), and with AlphaMissense or REVEL missing; each with Katz intervals.
It also counts, for every test set of benchmarks/results/acmg_calibration.json, the bands whose point
estimate and whose lower Katz bound reach the required ratio.

Writes benchmarks/results/calibration_by_inputs.{json,md}.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))

import calibrate_acmg as C
import train_model as T

OUT = Path(__file__).resolve().parent / "results"


def _f(x) -> float:
    """A ratio as a number ('inf' stays infinite; a missing one counts as not reaching)."""
    return float("-inf") if x is None else float("inf") if x == "∞" else float(x)


def main() -> int:
    rep: dict = {"patterns": {}, "bands_reaching": {}}
    j = json.loads((OUT / "acmg_calibration.json").read_text())
    cal = {k: {int(a): b for a, b in j["calibration"]["model"][k].items()} for k in ("pp3", "bp4")}
    for name, ev in j["evaluation"].items():
        bands = ev["bands"]
        rep["bands_reaching"][name] = {
            "bands": len(bands),
            "point_estimate": sum(_f(b["lr"]) >= b["required"] for b in bands.values()),
            "lower_bound": sum(_f(b["ci_katz"][0]) >= b["required"] for b in bands.values()),
        }
    c20 = C.pejaver_clinvar("2020")
    train_ids = set(pd.read_csv(T.OUT / "model_training_set.csv.gz", dtype={"clinvar_id": str})["clinvar_id"])
    vids = C._grch38_map_cached("2020")
    d = c20[[vids[int(i)][1] not in train_ids for i in c20["id"]]].copy()
    d["model"] = C.model_scores(d)
    have = {c: d[c].notna() for c in ("alphamissense", "revel", "cadd_phred", "popeve")}
    patterns = {
        "all four scores": have["alphamissense"] & have["revel"] & have["cadd_phred"] & have["popeve"],
        "popEVE missing (AlphaMissense and REVEL present)": have["alphamissense"] & have["revel"] & ~have["popeve"],
        "AlphaMissense or REVEL missing": ~(have["alphamissense"] & have["revel"]) & d["model"].notna(),
    }
    rng = np.random.default_rng(0)
    for name, m in patterns.items():
        g = d[m]
        pts = C.points(C._num_col(g, "model"), cal["pp3"], cal["bp4"])
        bands = C.band_lrs(pts, g["label"].to_numpy(), rng, n_boot=500)
        rep["patterns"][name] = {
            "pathogenic": int(g["label"].sum()),
            "benign": int((g["label"] == 0).sum()),
            "bands": {
                str(k): {x: b[x] for x in ("n_pathogenic", "n_benign", "lr", "ci_katz", "required")}
                for k, b in bands.items()
            },
        }
    OUT.mkdir(exist_ok=True)
    (OUT / "calibration_by_inputs.json").write_text(json.dumps(rep, indent=1, default=str))
    lines = ["# PP3/BP4 bands by missing predictor scores", "",
             "Produced by `benchmarks/calibration_by_inputs.py`.", "",
             "Bands whose point estimate / lower Katz bound reach the required ratio, by test set:", ""]  # fmt: skip
    for name, x in rep["bands_reaching"].items():
        lines.append(f"- {name}: {x['point_estimate']} / {x['lower_bound']} of {x['bands']}")
    for name, x in rep["patterns"].items():
        lines += ["", f"## ClinVar 2020 set, shipped model: {name} ({x['pathogenic']} P / {x['benign']} B)", "",
                  "| Band | P | B | LR (Katz 95%) | Required |", "|---|---|---|---|---|"]  # fmt: skip
        for k, b in x["bands"].items():
            ci = b["ci_katz"]
            lines.append(
                f"| {k} | {b['n_pathogenic']} | {b['n_benign']} | {b['lr']} ({ci[0]}–{ci[1]}) | {b['required']} |"
            )
    (OUT / "calibration_by_inputs.md").write_text("\n".join(lines) + "\n")
    print("\n".join(lines))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
