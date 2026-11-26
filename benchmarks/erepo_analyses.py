"""Secondary analyses of the comparison with ClinGen's expert panels (offline).

    python benchmarks/erepo_analyses.py

Reads benchmarks/results/erepo_validation_variants.csv.gz (written by validate_erepo.py) and reports:

1. Per panel: PP3/BP4 agreement (sensitivity and specificity), so that the pooled figures can be
   compared with the spread between panels.
2. Clustered intervals: a panel-level cluster bootstrap of the per-criterion agreement (panels apply
   their own specifications, so their classifications are not independent); for the AUC differences
   (model - REVEL, model - AlphaMissense) the variant-level DeLong test, a gene-level cluster bootstrap
   and the AUC within genes (pathogenic-benign pairs of the same gene only).
3. PP3/BP4 strength: exact agreement and linear-weighted kappa of the signed points (BP4 negative),
   and kappa of the direction alone, among the classifications whose panel evaluated PP3 or BP4.
4. Points (Tavtigian et al. 2020: pathogenic >= 10, likely pathogenic 6-9, VUS 0-5, likely benign -1 to
   -6, benign <= -7; BA1 alone gives benign): the panel's class rebuilt from its met codes; where that
   reproduces the panel's class, the class after replacing the panel's PP3/BP4 points by the tool's
   (an illustration of what the tool's strengths would do, not a recommendation); and the class the
   tool's automated codes would give alone.
5. "Not met" lists: some panels list every code they did not apply as not met; agreement without the
   classifications that list 15 or more, and BP4 specificity without those where the panel applied
   BA1 or BS1 (there BP4 is moot).
6. PM2: missed applications by gnomAD grpmax and sensitivity by approval period; PS1/PM5 applied by the
   tool, by the review stars and class of the ClinVar reference variant.
7. Incorporation bias: the panels' classes partly rest on PP3/BP4 (often from REVEL); AUCs on the
   variants whose rebuilt class does not change when PP3/BP4 are removed.

Writes benchmarks/results/erepo_analyses.{json,md}.
"""

from __future__ import annotations

import json
import re
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import chi2, norm
from sklearn.metrics import cohen_kappa_score, roc_auc_score

OUT = Path(__file__).resolve().parent / "results"
CLASSES = ["Benign", "Likely Benign", "Uncertain Significance", "Likely Pathogenic", "Pathogenic"]
PLP, BLB = {"Pathogenic", "Likely Pathogenic"}, {"Benign", "Likely Benign"}
N_BOOT = 10000
SCORES = ("model_score", "revel", "alphamissense")


def tier(points: float, ba1: bool) -> str:
    """ACMG class of a point total (Tavtigian et al. 2020); BA1 alone gives Benign."""
    if ba1:
        return "Benign"
    if points >= 10:
        return "Pathogenic"
    if points >= 6:
        return "Likely Pathogenic"
    if points >= 0:
        return "Uncertain Significance"
    return "Likely Benign" if points >= -6 else "Benign"


DEFAULT = {"PVS": 8, "PS": 4, "PM": 2, "PP": 1, "BS": 4, "BP": 1}
SUFFIX = {"Supporting": 1, "Moderate": 2, "Strong": 4, "Very": 8, "Very Strong": 8, "VeryStrong": 8}


def met_points(text: str) -> dict[str, int]:
    """The panel's met codes with their points: 'PVS1,PM2_Supporting,PP3_Moderate' -> {PVS1: 8, PM2: 1, PP3: 2}.

    A code without a suffix has its default ACMG strength (PVS 8, PS/BS 4, PM 2, PP/BP 1); BA1 is kept
    as stand-alone benign.
    """
    out: dict[str, int] = {}
    for c in (x.strip() for x in re.split(r"[,;]", text or "")):
        if not c:
            continue
        base, _, mod = c.partition("_")
        if base.startswith("BA"):
            out["BA1"] = 8
            continue
        kind = re.match(r"^[A-Z]+", base)
        pts = SUFFIX.get(mod.strip()) if mod else None
        out[base] = max(out.get(base, 0), pts or DEFAULT.get(kind.group(0) if kind else "", 0))
    return out


def signed_points(met: dict[str, int], drop: tuple[str, ...] = ()) -> tuple[float, bool]:
    """Sum of the met codes' points (benign codes negative) and whether BA1 is met."""
    total = 0.0
    for code, pts in met.items():
        if code in drop or code == "BA1":
            continue
        total += -pts if code.startswith("B") else pts
    return total, "BA1" in met


def rate(k: int, n: int) -> float | None:
    return round(k / n, 6) if n else None


def agreement(d: pd.DataFrame, code: str) -> dict:
    ev = d[d[f"exp_{code}"].notna()]
    exp, tool = ev[f"exp_{code}"] > 0, ev[f"tool_{code}"] > 0
    tp, fn, fp, tn = (int(x.sum()) for x in (exp & tool, exp & ~tool, ~exp & tool, ~exp & ~tool))
    return {"tp": tp, "fn": fn, "fp": fp, "tn": tn, "sensitivity": rate(tp, tp + fn), "specificity": rate(tn, tn + fp)}


def _pct_ci(vals: np.ndarray) -> list[float | None]:
    vals = vals[~np.isnan(vals)]
    if not len(vals):
        return [None, None]
    return [round(float(np.percentile(vals, 2.5)), 6), round(float(np.percentile(vals, 97.5)), 6)]


def cluster_agreement_ci(d: pd.DataFrame, by: str, code: str, rng: np.random.Generator) -> dict:
    """Sensitivity and specificity with 95% intervals from a bootstrap of whole clusters (`by`)."""
    ev = d[d[f"exp_{code}"].notna()]
    exp, tool = (ev[f"exp_{code}"] > 0).to_numpy(), (ev[f"tool_{code}"] > 0).to_numpy()
    cells = (
        pd.DataFrame(
            {"c": ev[by].to_numpy(), "tp": exp & tool, "fn": exp & ~tool, "fp": ~exp & tool, "tn": ~exp & ~tool}
        )
        .groupby("c")[["tp", "fn", "fp", "tn"]]
        .sum()
        .to_numpy(float)
    )
    pick = rng.integers(0, len(cells), (N_BOOT, len(cells)))
    tot = cells[pick].sum(axis=1)  # N_BOOT x 4
    with np.errstate(invalid="ignore", divide="ignore"):
        sens, spec = tot[:, 0] / (tot[:, 0] + tot[:, 1]), tot[:, 3] / (tot[:, 3] + tot[:, 2])
        ppv = tot[:, 0] / (tot[:, 0] + tot[:, 2])
    return {"sensitivity_ci": _pct_ci(sens), "specificity_ci": _pct_ci(spec), "ppv_ci": _pct_ci(ppv),
            "clusters": len(cells)}  # fmt: skip


def cluster_auc_diff_ci(d: pd.DataFrame, by: str, a: str, b: str, rng: np.random.Generator) -> list[float | None]:
    """95% interval of AUC(a) - AUC(b) from a bootstrap of whole clusters (`by`)."""
    y, sa, sb = d["y"].to_numpy(), d[a].to_numpy(float), d[b].to_numpy(float)
    groups = [np.asarray(ix) for ix in d.reset_index(drop=True).groupby(by).indices.values()]
    vals = np.full(N_BOOT, np.nan)
    for k in range(N_BOOT):
        i = np.concatenate([groups[j] for j in rng.integers(0, len(groups), len(groups))])
        if 0 < y[i].sum() < len(i):
            vals[k] = roc_auc_score(y[i], sa[i]) - roc_auc_score(y[i], sb[i])
    return _pct_ci(vals)


# --------------------------------------------------------------------------- AUC: DeLong and clusters


def _midrank(x: np.ndarray) -> np.ndarray:
    j = np.argsort(x)
    z = x[j]
    n = len(x)
    t = np.zeros(n)
    i = 0
    while i < n:
        k = i
        while k < n and z[k] == z[i]:
            k += 1
        t[i:k] = 0.5 * (i + k - 1) + 1
        i = k
    out = np.empty(n)
    out[j] = t
    return out


def delong(y: np.ndarray, a: np.ndarray, b: np.ndarray) -> dict:
    """Paired DeLong test of AUC(a) - AUC(b) (fast algorithm of Sun and Xu, 2014)."""
    pos, neg = y == 1, y == 0
    m, n = int(pos.sum()), int(neg.sum())
    aucs, v10, v01 = [], [], []
    for s in (a, b):
        tx, ty, tz = _midrank(s[pos]), _midrank(s[neg]), _midrank(np.concatenate([s[pos], s[neg]]))
        aucs.append(tz[:m].sum() / (m * n) - (m + 1) / (2 * n))
        v10.append((tz[:m] - tx) / n)
        v01.append(1 - (tz[m:] - ty) / m)
    s10, s01 = np.cov(np.array(v10)), np.cov(np.array(v01))
    var = (s10[0, 0] + s10[1, 1] - 2 * s10[0, 1]) / m + (s01[0, 0] + s01[1, 1] - 2 * s01[0, 1]) / n
    diff = aucs[0] - aucs[1]
    se = float(np.sqrt(var))
    return {
        "auc_a": round(float(aucs[0]), 6),
        "auc_b": round(float(aucs[1]), 6),
        "difference": round(float(diff), 6),
        "ci": [round(diff - 1.959964 * se, 6), round(diff + 1.959964 * se, 6)],
        "p": float(2 * norm.sf(abs(diff) / se)) if se > 0 else None,
    }


def within_gene_diff_ci(d: pd.DataFrame, a: str, b: str, rng: np.random.Generator) -> list[float | None]:
    """95% interval of within-gene AUC(a) - AUC(b) from a bootstrap of whole genes."""
    per = []
    for _, g in d.groupby("gene"):
        p, q = g[g["y"] == 1], g[g["y"] == 0]
        if len(p) and len(q):
            cells = []
            for s in (a, b):
                diff = p[s].to_numpy(float)[:, None] - q[s].to_numpy(float)[None, :]
                cells.append((diff > 0).sum() + 0.5 * (diff == 0).sum())
            per.append((cells[0], cells[1], len(p) * len(q)))
    arr = np.array(per, float)
    pick = rng.integers(0, len(arr), (N_BOOT, len(arr)))
    tot = arr[pick].sum(axis=1)
    return _pct_ci((tot[:, 0] - tot[:, 1]) / tot[:, 2])


def within_gene_auc(d: pd.DataFrame, score: str) -> float | None:
    """AUC over pathogenic-benign pairs of the same gene only (ties count half)."""
    conc = pairs = 0.0
    for _, g in d.groupby("gene"):
        p, b = g.loc[g["y"] == 1, score].to_numpy(), g.loc[g["y"] == 0, score].to_numpy()
        if len(p) and len(b):
            diff = p[:, None] - b[None, :]
            conc += (diff > 0).sum() + 0.5 * (diff == 0).sum()
            pairs += diff.size
    return round(conc / pairs, 6) if pairs else None


def auc_comparisons(lab: pd.DataFrame, rng: np.random.Generator) -> dict:
    d = lab.dropna(subset=list(SCORES))
    y = d["y"].to_numpy()
    out = {
        "n": len(d),
        "pathogenic": int(y.sum()),
        "benign": int((y == 0).sum()),
        "genes": int(d["gene"].nunique()),
        "auc": {s: round(float(roc_auc_score(y, d[s])), 6) for s in SCORES},
        "within_gene_auc": {s: within_gene_auc(d, s) for s in SCORES},
        "within_gene_pairs": int(sum(((g["y"] == 1).sum() * (g["y"] == 0).sum()) for _, g in d.groupby("gene"))),
        "within_gene_model_minus_revel_ci": within_gene_diff_ci(d, "model_score", "revel", rng),
    }
    for other in ("revel", "alphamissense"):
        dl = delong(y, d["model_score"].to_numpy(float), d[other].to_numpy(float))
        out[f"model - {other}"] = {
            **dl,
            "gene_cluster_ci": cluster_auc_diff_ci(d, "gene", "model_score", other, rng),
            "panel_cluster_ci": cluster_auc_diff_ci(d, "panel", "model_score", other, rng),
        }
    return out


# --------------------------------------------------------------------------- main


def main(argv: list[str] | None = None) -> int:
    import argparse

    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--table", type=Path, default=OUT / "erepo_validation_variants.csv.gz")
    ap.add_argument("--out", type=Path, default=OUT)
    args = ap.parse_args(argv)
    df = pd.read_csv(args.table, dtype={"clinvar_id": str})
    df["met"] = df["panel_codes"].fillna("").map(met_points)
    df["n_not_met"] = df["panel_not_met"].fillna("").map(lambda s: len([x for x in s.split(",") if x]))
    rng = np.random.default_rng(0)
    rep: dict = {"classifications": len(df), "panels": int(df["panel"].nunique()), "genes": int(df["gene"].nunique())}

    # 1-2. per panel and clustered agreement
    per_panel = {}
    for panel, g in df.groupby("panel"):
        per_panel[panel] = {
            "classifications": len(g),
            **{c: agreement(g, c) for c in ("PP3", "BP4")},
        }
    rep["per_panel"] = per_panel
    rep["agreement_panel_cluster_ci"] = {}
    for c in ("PP3", "BP4", "PM2", "PS1", "PM5", "BA1"):
        a = agreement(df, c)
        rep["agreement_panel_cluster_ci"][c] = {**a, **cluster_agreement_ci(df, "panel", c, rng)}

    lab = df[df["assertion"].isin(PLP | BLB)].copy()
    lab["y"] = lab["assertion"].isin(PLP).astype(int)
    rep["auc"] = {
        "all": auc_comparisons(lab, rng),
        "gene never seen": auc_comparisons(lab[lab["overlap"] == "gene never seen"], rng),
    }

    # 3. strength agreement
    ev = df[df["exp_PP3"].notna() | df["exp_BP4"].notna()].copy()
    panel_s = np.where(ev["exp_PP3"].fillna(0) > 0, ev["exp_PP3"].fillna(0), -ev["exp_BP4"].fillna(0)).astype(int)
    tool_s = np.where(ev["tool_PP3"] > 0, ev["tool_PP3"], -ev["tool_BP4"]).astype(int)
    labels = sorted(set(panel_s) | set(tool_s))
    rep["strength"] = {
        "classifications": len(ev),
        "exact_agreement": rate(int((panel_s == tool_s).sum()), len(ev)),
        "weighted_kappa_linear": round(float(cohen_kappa_score(panel_s, tool_s, labels=labels, weights="linear")), 4),
        "kappa_direction": round(float(cohen_kappa_score(np.sign(panel_s), np.sign(tool_s))), 4),
        "pp3_panel_met": {
            "n": int((ev["exp_PP3"] > 0).sum()),
            "same": int(((ev["exp_PP3"] > 0) & (ev["exp_PP3"] == ev["tool_PP3"])).sum()),
            "higher": int(((ev["exp_PP3"] > 0) & (ev["tool_PP3"] > ev["exp_PP3"])).sum()),
            "lower_not_none": int(
                ((ev["exp_PP3"] > 0) & (ev["tool_PP3"] > 0) & (ev["tool_PP3"] < ev["exp_PP3"])).sum()
            ),
            "none": int(((ev["exp_PP3"] > 0) & (ev["tool_PP3"] == 0)).sum()),
        },
        "bp4_panel_met": {
            "n": int((ev["exp_BP4"] > 0).sum()),
            "same": int(((ev["exp_BP4"] > 0) & (ev["exp_BP4"] == ev["tool_BP4"])).sum()),
            "higher": int(((ev["exp_BP4"] > 0) & (ev["tool_BP4"] > ev["exp_BP4"])).sum()),
            "lower_not_none": int(
                ((ev["exp_BP4"] > 0) & (ev["tool_BP4"] > 0) & (ev["tool_BP4"] < ev["exp_BP4"])).sum()
            ),
            "none": int(((ev["exp_BP4"] > 0) & (ev["tool_BP4"] == 0)).sum()),
        },
        "panel_pp3_strength": {int(k): int(v) for k, v in ev.loc[ev["exp_PP3"] > 0, "exp_PP3"].value_counts().items()},
        "panels_applying_pp3_above_supporting": int(ev.loc[ev["exp_PP3"] > 1, "panel"].nunique()),
    }

    # 4. points
    rebuilt = []
    for _, r in df.iterrows():
        met = r["met"]
        base, ba1 = signed_points(met, drop=("PP3", "BP4"))
        full, _ = signed_points(met)
        tool_pp = (r["tool_PP3"] or 0) - (r["tool_BP4"] or 0)
        auto = tool_pp + (r["tool_PS1"] or 0) + (r["tool_PM5"] or 0) + (r["tool_PM2"] or 0)
        rebuilt.append(
            {
                "rebuilt": tier(full, ba1),
                "without_pp3_bp4": tier(base, ba1),
                "with_tool_pp3_bp4": tier(base + tool_pp, ba1),
                "automated_only": tier(auto, bool(r["tool_BA1"])),
            }
        )
    pts = pd.concat([df.reset_index(drop=True), pd.DataFrame(rebuilt)], axis=1)
    ok = pts[pts["rebuilt"] == pts["assertion"]]
    rep["points"] = {
        "rebuilt_matches_panel": len(ok),
        "of": len(pts),
        "substitution": {
            f"{a} -> {b}": int(n)
            for (a, b), n in ok.groupby(["assertion", "with_tool_pp3_bp4"]).size().items()
            if a != b
        },
        "substitution_by_panel_class": {a: int(n) for a, n in ok["assertion"].value_counts().items()},
        "automated_only": {
            f"{a} | {b}": int(n) for (a, b), n in pts.groupby(["assertion", "automated_only"]).size().items()
        },
    }

    # 5. "not met" lists
    blanket = pts["n_not_met"] >= 15
    rep["not_met_lists"] = {
        "classifications_with_15_or_more": int(blanket.sum()),
        "by_panel_share": {
            p: round(float(g.mean()), 3) for p, g in blanket.groupby(pts["panel"]) if g.mean() >= 0.5 and len(g) >= 20
        },
        "agreement_without_them": {c: agreement(pts[~blanket], c) for c in ("PP3", "BP4", "PM2", "PS1", "PM5", "BA1")},
        "bp4_specificity_without_ba1_bs1": agreement(pts[~pts["met"].map(lambda m: "BA1" in m or "BS1" in m)], "BP4")[
            "specificity"
        ],
    }

    # 6. PM2 and PS1/PM5 details
    miss = pts[(pts["exp_PM2"] > 0) & (pts["tool_PM2"] == 0)]
    gm = miss["grpmax"]
    year = pts["approved"].str[:4].astype(int)
    rep["pm2"] = {
        "missed": len(miss),
        "missed_by_grpmax": {
            "no frequency": int(gm.isna().sum()),
            "< 1e-5": int((gm < 1e-5).sum()),
            "1e-5 to 1e-4": int(((gm >= 1e-5) & (gm < 1e-4)).sum()),
            ">= 1e-4": int((gm >= 1e-4).sum()),
        },
        "missed_grpmax_median": None if gm.dropna().empty else float(gm.median()),
        "sensitivity_by_period": {
            name: agreement(pts[m], "PM2")["sensitivity"]
            for name, m in (("2021 or earlier", year <= 2021), ("2022-2023", (year >= 2022) & (year <= 2023)),
                            ("2024 or later", year >= 2024))
        },
    }  # fmt: skip
    rep["pm2"]["within_panel_2024_vs_2021"] = pm2_period_mh(pts)
    ps = pts[(pts["tool_PS1"] > 0) | (pts["tool_PM5"] > 0)].copy()
    # the first reference listed is the strongest (residue.py sorts by strength): its level and stars
    first = (
        ps["tool_ps1_pm5"]
        .astype(str)
        .str.extract(r"\((?:same change by )?[^:]+: (Pathogenic|Likely pathogenic), (\d)★")
    )
    ps["stars"] = first[1]
    ps["ref_class"] = np.where(first[0] == "Likely pathogenic", "LP", "P")
    out = {}
    for c in ("PS1", "PM5"):
        d = ps[(ps[f"tool_{c}"] > 0) & ps[f"exp_{c}"].notna()]
        out[c] = {
            f"{s} stars, {k}": {
                "panel met": int((g[f"exp_{c}"] > 0).sum()),
                "panel not met": int((g[f"exp_{c}"] == 0).sum()),
            }
            for (s, k), g in d.groupby(["stars", "ref_class"])
        }
    rep["ps1_pm5_by_reference"] = out

    # 7. incorporation bias
    stable = pts[(pts["rebuilt"] == pts["assertion"]) & (pts["without_pp3_bp4"] == pts["assertion"])]
    sl = stable[stable["assertion"].isin(PLP | BLB)].copy()
    sl["y"] = sl["assertion"].isin(PLP).astype(int)
    rep["incorporation_bias"] = auc_comparisons(sl, rng)

    # SpliceAI and BP4 (v2.5 withholds BP4 above 0.1)
    if "spliceai_max" in pts:
        rep["bp4_with_spliceai_above_0.1"] = int(((pts["tool_BP4"] > 0) & (pts["spliceai_max"] > 0.1)).sum())

    args.out.mkdir(exist_ok=True)
    (args.out / "erepo_analyses.json").write_text(json.dumps(rep, indent=1, default=str))
    (args.out / "erepo_analyses.md").write_text(markdown(rep))
    print(markdown(rep))
    return 0


def pm2_period_mh(pts: pd.DataFrame) -> dict:
    """Mantel-Haenszel odds ratio of the tool's PM2 among the panels' PM2, from 2024 against up to 2021.

    Stratified by panel: does the rise of PM2 agreement with approval date hold within panels?
    """
    d = pts[pts["exp_PM2"] > 0].copy()
    late, early = d["approved"] >= "2024-01-01", d["approved"] <= "2021-12-31"
    num = den = a_sum = e_sum = var_sum = 0.0
    rbg = [0.0, 0.0, 0.0]  # Robins-Breslow-Greenland variance terms
    k = 0
    for _, g in d[late | early].groupby("panel"):
        lt, hit = g["approved"] >= "2024-01-01", g["tool_PM2"] > 0
        a, b, c, e = int((lt & hit).sum()), int((lt & ~hit).sum()), int((~lt & hit).sum()), int((~lt & ~hit).sum())
        n = a + b + c + e
        if min(a + b, c + e, a + c, b + e) == 0:
            continue
        k += 1
        num += a * e / n
        den += b * c / n
        pp, qq = (a + e) / n, (b + c) / n
        rbg[0] += pp * a * e / n
        rbg[1] += pp * b * c / n + qq * a * e / n
        rbg[2] += qq * b * c / n
        a_sum += a
        e_sum += (a + b) * (a + c) / n
        var_sum += (a + b) * (c + e) * (a + c) * (b + e) / (n * n * (n - 1))
    stat = (abs(a_sum - e_sum) - 0.5) ** 2 / var_sum
    or_mh = num / den
    se = np.sqrt(rbg[0] / (2 * num**2) + rbg[1] / (2 * num * den) + rbg[2] / (2 * den**2))
    ci = [round(float(np.exp(np.log(or_mh) + z * se)), 3) for z in (-1.96, 1.96)]
    return {"panels": k, "odds_ratio": round(or_mh, 3), "ci_rbg": ci, "p": float(chi2.sf(stat, 1))}


def _p(x: float | None) -> str:
    return "—" if x is None else f"{100 * x:.1f}%"


def markdown(rep: dict) -> str:
    L = [
        "# eRepo comparison: secondary analyses",
        "",
        f"{rep['classifications']:,} classifications, {rep['panels']} panels, {rep['genes']} genes. Produced by "
        "`benchmarks/erepo_analyses.py` from `erepo_validation_variants.csv.gz`.",
        "",
        "## Agreement with panel-clustered 95% intervals",
        "",
        "| Criterion | Positive agreement (panel-cluster CI) | Negative agreement (panel-cluster CI) |",
        "|---|---|---|",
    ]
    for c, a in rep["agreement_panel_cluster_ci"].items():
        s, sp = a["sensitivity_ci"], a["specificity_ci"]
        L.append(
            f"| {c} | {_p(a['sensitivity'])} ({_p(s[0])}–{_p(s[1])}) | "
            f"{_p(a['specificity'])} ({_p(sp[0])}–{_p(sp[1])}) |"
        )
    for name, a in rep["auc"].items():
        L += ["", f"## AUC against the panel's class: {name}", ""]
        L.append(f"n = {a['n']} ({a['pathogenic']} P/LP, {a['benign']} B/LB), {a['genes']} genes.")
        L.append(f"- AUC: {json.dumps(a['auc'])}; within genes ({a['within_gene_pairs']:,} pairs): "
                 f"{json.dumps(a['within_gene_auc'])}")  # fmt: skip
        for k in ("model - revel", "model - alphamissense"):
            x = a[k]
            L.append(
                f"- {k}: {x['difference']:+.4f} (DeLong {x['ci'][0]:+.4f} to {x['ci'][1]:+.4f}, p = {x['p']:.2g}; "
                f"gene-cluster bootstrap {x['gene_cluster_ci'][0]:+.4f} to {x['gene_cluster_ci'][1]:+.4f})"
            )
    s = rep["strength"]
    L += [
        "",
        "## PP3/BP4 strength",
        "",
        f"- Classifications with PP3 or BP4 evaluated: {s['classifications']:,}; exact agreement of signed points "
        f"{_p(s['exact_agreement'])}; linear-weighted kappa {s['weighted_kappa_linear']}; kappa of direction "
        f"{s['kappa_direction']}.",
        f"- Panel applied PP3: {json.dumps(s['pp3_panel_met'])}; strengths {json.dumps(s['panel_pp3_strength'])}; "
        f"panels applying PP3 above supporting: {s['panels_applying_pp3_above_supporting']}.",
        f"- Panel applied BP4: {json.dumps(s['bp4_panel_met'])}.",
    ]
    p = rep["points"]
    L += [
        "",
        "## Points",
        "",
        f"- Panel class rebuilt from its met codes: {p['rebuilt_matches_panel']:,} of {p['of']:,} reproduce it.",
        f"- After replacing the panel's PP3/BP4 by the tool's (among those): {json.dumps(p['substitution'])}; "
        f"by panel class: {json.dumps(p['substitution_by_panel_class'])}.",
        f"- Class from the tool's automated codes alone (panel | automated): {json.dumps(p['automated_only'])}.",
    ]
    n = rep["not_met_lists"]
    L += [
        "",
        "## 'Not met' lists",
        "",
        f"- Classifications listing >= 15 codes as not met: {n['classifications_with_15_or_more']:,}; panels where "
        f"they are the majority: {json.dumps(n['by_panel_share'])}.",
        "- Agreement without them: "
        + "; ".join(
            f"{c} {_p(a['sensitivity'])} / {_p(a['specificity'])}" for c, a in n["agreement_without_them"].items()
        ),
        f"- BP4 specificity without classifications where the panel applied BA1 or BS1: "
        f"{_p(n['bp4_specificity_without_ba1_bs1'])}.",
        "",
        "## PM2 and PS1/PM5",
        "",
        f"- PM2 missed: {rep['pm2']['missed']:,}; by grpmax {json.dumps(rep['pm2']['missed_by_grpmax'])}; median "
        f"grpmax {rep['pm2']['missed_grpmax_median']}; sensitivity by period "
        f"{json.dumps(rep['pm2']['sensitivity_by_period'])}.",
        "- PS1/PM5 applied by the tool, by the reference's stars and class: "
        f"{json.dumps(rep['ps1_pm5_by_reference'])}.",
    ]
    ib = rep["incorporation_bias"]
    L += [
        "",
        "## Incorporation bias",
        "",
        f"Variants whose rebuilt class does not depend on PP3/BP4: {ib['n']} ({ib['pathogenic']} P/LP, "
        f"{ib['benign']} B/LB). AUC {json.dumps(ib['auc'])}; model - REVEL {ib['model - revel']['difference']:+.4f} "
        f"(DeLong {ib['model - revel']['ci'][0]:+.4f} to {ib['model - revel']['ci'][1]:+.4f}; gene-cluster "
        f"{ib['model - revel']['gene_cluster_ci'][0]:+.4f} to {ib['model - revel']['gene_cluster_ci'][1]:+.4f}).",
    ]
    if "bp4_with_spliceai_above_0.1" in rep:
        L.append(f"\nBP4 given with SpliceAI above 0.1: {rep['bp4_with_spliceai_above_0.1']}.")
    L.append("")
    return "\n".join(L)


if __name__ == "__main__":
    raise SystemExit(main())
