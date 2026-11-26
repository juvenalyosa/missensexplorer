"""Are the panels' decisions reproduced by their own specifications? A consistency check (offline).

    python benchmarks/erepo_vcep_rules.py

validate_erepo.py compares the tool's generic criteria with ClinGen expert panels' decisions. Panels apply
their own specifications (ClinGen CSpec). Here the PP3/BP4, PM2 and BA1 rules of the ten panels with most
classifications in that comparison (src/missensexplorer/vcep_rules.py, encoded by hand from CSpec) are
applied to the same variants and compared with the panels' decisions, next to the tool's generic criteria
and next to the tool's PP3/BP4 capped at the highest strength each specification allows (the structured
part of CSpec, src/missensexplorer/vcep_params.py).

This is a consistency check, not a test of which strength is right: panels wrote these rules and then
applied them, so high agreement shows that the disagreements with the tool are systematic. Limits, stated
in the report:
- the ten panels were chosen by size, and most have simple single-threshold REVEL rules;
- the rules are those of the current specification versions, applied to decisions made since 2016, so
  agreement is also given by approval period;
- frequency thresholds that panels apply to the filtering allele frequency (FAF) are applied to gnomAD's
  grpmax allele frequency, which is at least the FAF: PM2 is given less often and BA1 more often than
  with the FAF, so BA1 positive agreement is an upper bound;
- homozygote counts, coverage requirements, splicing tools other than SpliceAI and exception lists are
  not applied.

Intervals: 95% percentile bootstrap over panels (10 clusters, 2,000 resamples).

Reads benchmarks/results/erepo_validation_variants.csv.gz; writes benchmarks/results/erepo_vcep_rules.{json,md}
and the per-classification codes (benchmarks/results/erepo_vcep_rules_variants.csv.gz).
"""

# ruff: noqa: E501  (report tables are kept on one line each)
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from missensexplorer import vcep
from missensexplorer.vcep_rules import RULES, panel_rule_codes

OUT = Path(__file__).resolve().parent / "results"
CRITERIA = ("PP3", "BP4", "PM2", "BA1")
PERIODS = (
    ("2016-2021", None, "2021-12-31"),
    ("2022-2023", "2022-01-01", "2023-12-31"),
    ("2024-2026", "2024-01-01", None),
)
B = 10000


def agree(d: pd.DataFrame, exp: str, got: str) -> dict:
    """Positive and negative agreement of `got` with the panels' `exp` (evaluated classifications only)."""
    ev = d[d[exp].notna()]
    e, g = ev[exp] > 0, ev[got] > 0
    tp, fn, fp, tn = (int(x.sum()) for x in (e & g, e & ~g, ~e & g, ~e & ~g))
    r = lambda a, b: round(a / b, 6) if b else None  # noqa: E731
    return {"evaluated": len(ev), "tp": tp, "fn": fn, "fp": fp, "tn": tn,
            "positive_agreement": r(tp, tp + fn), "negative_agreement": r(tn, tn + fp)}  # fmt: skip


def signed(d: pd.DataFrame, prefix: str) -> pd.Series:
    """PP3 points as positive, BP4 points as negative (0 when neither)."""
    pp3, bp4 = d[f"{prefix}PP3"].fillna(0), d[f"{prefix}BP4"].fillna(0)
    return pd.Series(np.where(pp3 > 0, pp3, -bp4), index=d.index)


def cap(row: dict) -> tuple[int, int]:
    """Highest PP3 and BP4 points the gene's specification allows for a missense variant (CSpec)."""
    best = {"PP3": 0, "BP4": 0}
    for s in vcep.specifications(row["gene"]):
        if not str(s["panel"]).startswith(RULES[row["gene"]]["panel"]):
            continue
        for c in best:
            if c in s.get("not_missense", {}):
                continue
            pts = [vcep.POINTS[a] for a in s["allowed"].get(c) or []]
            best[c] = max([best[c], *pts])
    return best["PP3"], best["BP4"]


def strength(d: pd.DataFrame) -> dict:
    """Exact agreement of the PP3/BP4 strength with the panels', for the tool, the capped tool and the rules."""
    panel = signed(d, "exp_")
    out = {}
    subsets = {
        # the metric of the whole expert-panel comparison: PP3 strength where the panel applied PP3
        "panel_pp3_met": d["exp_PP3"].fillna(0) > 0,
        "panel_pp3_or_bp4_met": (d["exp_PP3"].fillna(0) > 0) | (d["exp_BP4"].fillna(0) > 0),
        "pp3_or_bp4_evaluated": d["exp_PP3"].notna() | d["exp_BP4"].notna(),
    }
    for name, m in subsets.items():
        out[name] = {"n": int(m.sum())}
        for who in ("tool", "tool_capped", "rule"):
            out[name][who] = round(float((signed(d[m], f"{who}_") == panel[m]).mean()), 6) if m.any() else None
    return out


PICKS = np.random.default_rng(1).integers(0, 10, (B, 10))  # the same panel resamples for every statistic


def _num(x: object) -> float | None:
    return None if x is None or (isinstance(x, float) and np.isnan(x)) else float(x)  # type: ignore[arg-type]


def tool_row(r: dict) -> dict:
    """The fields of a results row that vcep.notes() reads, rebuilt from the eRepo per-variant table."""
    pts = int(r["tool_PP3"]) if r["tool_PP3"] > 0 else -int(r["tool_BP4"]) if r["tool_BP4"] > 0 else 0
    return {"gene": r["gene"], "acmg_points": pts, "acmg_ps1_pm5": r["tool_ps1_pm5"], "acmg_pm2": r["tool_pm2"],
            "gnomad_popmax_af": _num(r["grpmax"]), "population_frequency": r["population_frequency"],
            "revel": _num(r["revel"]), "spliceai_max": _num(r["spliceai_max"])}  # fmt: skip


def notes_on_erepo(df: pd.DataFrame, d: pd.DataFrame) -> dict:
    """What the tool's expert-panel notes (vcep.notes, as predict.py calls it) say on the eRepo classifications."""
    notes = {i: " | ".join(vcep.notes(tool_row(r))) for i, r in df.to_dict("index").items()}
    text = pd.Series(notes)
    # the ten panels: the note gives the panel rules' PP3/BP4, which differs from the tool's here
    differ = signed(d, "rule_") != signed(d, "tool_")
    given = text[d.index].str.contains("its rules give")
    ba1 = (df["exp_BA1"] > 0) & (df["tool_BA1"] == 0)  # the panel applied BA1 and the tool did not
    cites = text.str.contains(r"BA1 threshold|BA1 rule|BA1 only", regex=True)
    ten = df.index.isin(d.index)
    return {
        "ten_panels_rules_pp3_bp4_differ_from_tool": [int((differ & given).sum()), len(d)],
        "ba1_panel_met_tool_not_note_cites_panel_ba1": {
            "all_panels": [int((ba1 & cites).sum()), int(ba1.sum())],
            "ten_panels": [int((ba1 & cites & ten).sum()), int((ba1 & ten).sum())],
        },
        "classifications_with_a_note": [int(text.ne("").sum()), len(df)],
    }


def bootstrap(d: pd.DataFrame, num, den) -> list[float] | None:
    """95% percentile interval of sum(num)/sum(den) over panels resampled with replacement.

    `num` and `den` map one panel's rows to a count, so each resample is a sum of per-panel counts.
    """
    groups = [g for _, g in d.groupby("panel")]
    assert len(groups) == PICKS.shape[1]
    n = np.array([num(g) for g in groups], float)
    m = np.array([den(g) for g in groups], float)
    tot_n, tot_m = n[PICKS].sum(1), m[PICKS].sum(1)
    ok = tot_m > 0
    vals = tot_n[ok] / tot_m[ok]
    return [round(float(x), 4) for x in np.percentile(vals, [2.5, 97.5])] if len(vals) else None


def _agree_counts(c: str, who: str, k: str):
    """(numerator, denominator) functions of positive or negative agreement, for bootstrap()."""
    e, w = f"exp_{c}", f"{who}_{c}"
    if k == "positive_agreement":
        return (lambda g: int(((g[e] > 0) & (g[w] > 0)).sum()), lambda g: int((g[e] > 0).sum()))
    return (lambda g: int(((g[e] == 0) & (g[w] == 0)).sum()), lambda g: int((g[e] == 0).sum()))


def _subset(g: pd.DataFrame, name: str) -> pd.Series:
    if name == "panel_pp3_met":
        return g["exp_PP3"].fillna(0) > 0
    if name == "panel_pp3_or_bp4_met":
        return (g["exp_PP3"].fillna(0) > 0) | (g["exp_BP4"].fillna(0) > 0)
    return g["exp_PP3"].notna() | g["exp_BP4"].notna()


def _strength_counts(name: str, who: str):
    """(numerator, denominator) functions of exact strength agreement in a subset, for bootstrap()."""

    def num(g: pd.DataFrame) -> int:
        x = g[_subset(g, name)]
        return int((signed(x, f"{who}_") == signed(x, "exp_")).sum())

    return num, lambda g: int(_subset(g, name).sum())


def main() -> int:
    df = pd.read_csv(OUT / "erepo_validation_variants.csv.gz", dtype={"clinvar_id": str})
    # a rule applies to its own panel's classifications only (RYR1 is also curated by other panels)
    mask = [g in RULES and str(p).startswith(RULES[g]["panel"]) for g, p in zip(df["gene"], df["panel"], strict=True)]
    d = df[mask].copy()
    d["n_not_met"] = d["panel_not_met"].fillna("").map(lambda x: len([c for c in str(x).split(",") if c]))
    codes = [panel_rule_codes(r, RULES[r["gene"]], af_key="grpmax", pm2_key="tool_pm2") for r in d.to_dict("records")]
    for c in CRITERIA:
        d[f"rule_{c}"] = [x[c] for x in codes]
    caps = [cap(r) for r in d.to_dict("records")]
    d["tool_capped_PP3"] = np.minimum(d["tool_PP3"].fillna(0), [c[0] for c in caps])
    d["tool_capped_BP4"] = np.minimum(d["tool_BP4"].fillna(0), [c[1] for c in caps])
    rep: dict = {"classifications": len(d), "of": len(df), "genes": int(d["gene"].nunique()),
                 "panels": sorted(d["panel"].unique()), "criteria": {}, "intervals": {}}  # fmt: skip
    for c in CRITERIA:
        rep["criteria"][c] = {
            "tool": agree(d, f"exp_{c}", f"tool_{c}"),
            "panel_rules": agree(d, f"exp_{c}", f"rule_{c}"),
        }
        for who in ("tool", "rule"):
            for k in ("positive_agreement", "negative_agreement"):
                rep["intervals"][f"{c}_{who}_{k}"] = bootstrap(d, *_agree_counts(c, who, k))
    rep["strength_exact_agreement"] = strength(d)
    for name in ("panel_pp3_met", "pp3_or_bp4_evaluated"):
        for who in ("tool", "tool_capped", "rule"):
            rep["intervals"][f"strength_{name}_{who}"] = bootstrap(d, *_strength_counts(name, who))
    myeloid = d["panel"].str.startswith("Myeloid")
    against = (d["exp_PM2"] == 0) & (d["rule_PM2"] > 0)
    rep["pm2_without_myeloid"] = {
        "tool": agree(d[~myeloid], "exp_PM2", "tool_PM2"),
        "panel_rules": agree(d[~myeloid], "exp_PM2", "rule_PM2"),
        "rule_conflicts_total": int(against.sum()),
        "rule_conflicts_myeloid": int((against & myeloid).sum()),
        "rule_conflicts_with_blanket_list": int((against & (d["n_not_met"] >= 15)).sum()),
    }
    year = d["approved"].str[:4]
    rep["pm2_myeloid_rule_conflicts_by_year"] = {
        y: {"conflicts": int((against & myeloid & (year == y)).sum()),
            "panel_not_met": int(((d["exp_PM2"] == 0) & myeloid & (year == y)).sum())}
        for y in sorted(year[myeloid].unique())
    }  # fmt: skip
    ba1_against = d[(d["exp_BA1"] == 0) & (d["rule_BA1"] > 0)]
    rep["ba1_rule_against_panel"] = {
        "n": len(ba1_against),
        "panel_used_bs1": int(ba1_against["panel_codes"].fillna("").str.contains(r"\bBS1").sum()),
    }
    rep["by_period"] = {}
    rep["by_period_without_phenylketonuria"] = {}
    pku = d["panel"].str.startswith("Phenylketonuria")
    for name, lo, hi in PERIODS:
        m = pd.Series(True, index=d.index)
        if lo:
            m &= d["approved"] >= lo
        if hi:
            m &= d["approved"] <= hi
        g = d[m]
        rep["by_period"][name] = {
            "classifications": len(g),
            "PP3": {"tool": agree(g, "exp_PP3", "tool_PP3"), "panel_rules": agree(g, "exp_PP3", "rule_PP3")},
            "BP4": {"tool": agree(g, "exp_BP4", "tool_BP4"), "panel_rules": agree(g, "exp_BP4", "rule_BP4")},
            "strength": strength(g)["pp3_or_bp4_evaluated"],
            "strength_panel_pp3_met": strength(g)["panel_pp3_met"],
        }
        h = d[m & ~pku]
        rep["by_period_without_phenylketonuria"][name] = {
            "classifications": len(h),
            "PP3": {"tool": agree(h, "exp_PP3", "tool_PP3"), "panel_rules": agree(h, "exp_PP3", "rule_PP3")},
            "strength": strength(h)["pp3_or_bp4_evaluated"],
            "strength_panel_pp3_met": strength(h)["panel_pp3_met"],
        }
    rep["by_panel"] = {
        p: {"classifications": len(g),
            **{c: {"tool": agree(g, f"exp_{c}", f"tool_{c}"), "panel_rules": agree(g, f"exp_{c}", f"rule_{c}")} for c in CRITERIA},
            "strength": strength(g)["pp3_or_bp4_evaluated"]}
        for p, g in d.groupby("panel")
    }  # fmt: skip
    rep["notes_on_erepo"] = notes_on_erepo(df, d)
    OUT.mkdir(exist_ok=True)
    keep = ["clinvar_id", "gene", "panel", "approved", "n_not_met",
            *[f"{w}_{c}" for c in CRITERIA for w in ("exp", "tool", "rule")], "tool_capped_PP3", "tool_capped_BP4"]  # fmt: skip
    d[keep].to_csv(OUT / "erepo_vcep_rules_variants.csv.gz", index=False, compression="gzip")
    (OUT / "erepo_vcep_rules.json").write_text(json.dumps(rep, indent=1))
    lines = report(rep)
    (OUT / "erepo_vcep_rules.md").write_text("\n".join(lines))
    print("\n".join(lines[:24]))
    return 0


def report(rep: dict) -> list[str]:
    pct = lambda x: "—" if x is None else f"{100 * x:.1f}%"  # noqa: E731

    def ci(k: str) -> str:
        v = rep["intervals"].get(k)
        return "" if not v else f" ({100 * v[0]:.1f}–{100 * v[1]:.1f})"

    lines = ["# The ten largest panels' own rules against their decisions (a consistency check)", "",
             "Produced by `benchmarks/erepo_vcep_rules.py` (rules in `src/missensexplorer/vcep_rules.py`, encoded by "
             "hand from ClinGen CSpec; see the script for the limits). Intervals: bootstrap over the 10 panels.", "",
             "| Criterion | Evaluated | Tool: positive / negative agreement | Panel rules: positive / negative agreement |",
             "|---|---|---|---|"]  # fmt: skip
    for c, x in rep["criteria"].items():
        t, r = x["tool"], x["panel_rules"]
        lines.append(f"| {c} | {t['evaluated']:,} | {pct(t['positive_agreement'])}{ci(f'{c}_tool_positive_agreement')} / "
                     f"{pct(t['negative_agreement'])}{ci(f'{c}_tool_negative_agreement')} | "
                     f"{pct(r['positive_agreement'])}{ci(f'{c}_rule_positive_agreement')} / "
                     f"{pct(r['negative_agreement'])}{ci(f'{c}_rule_negative_agreement')} |")  # fmt: skip
    lines += ["", "Exact agreement of the PP3/BP4 strength with the panels':", "",
              "| Classifications | n | Tool | Tool capped at the specification's maximum | Panel rules |",
              "|---|---|---|---|---|"]  # fmt: skip
    for name, x in rep["strength_exact_agreement"].items():
        lines.append(f"| {name} | {x['n']:,} | {pct(x['tool'])}{ci(f'strength_{name}_tool')} | "
                     f"{pct(x['tool_capped'])}{ci(f'strength_{name}_tool_capped')} | "
                     f"{pct(x['rule'])}{ci(f'strength_{name}_rule')} |")  # fmt: skip
    w = rep["pm2_without_myeloid"]
    lines += ["", f"PM2 without the Myeloid Malignancy panel: panel rules {pct(w['panel_rules']['positive_agreement'])} / "
              f"{pct(w['panel_rules']['negative_agreement'])}; tool {pct(w['tool']['positive_agreement'])} / "
              f"{pct(w['tool']['negative_agreement'])}. Rule PM2 against the panel: {w['rule_conflicts_total']} "
              f"({w['rule_conflicts_myeloid']} Myeloid Malignancy; {w['rule_conflicts_with_blanket_list']} in "
              "classifications listing 15 or more criteria as not met).", "", "By approval period:", "", "```",
              json.dumps(rep["by_period"], indent=1), "```", "", "By panel:", "", "```",
              json.dumps(rep["by_panel"], indent=1), "```", ""]  # fmt: skip
    return lines


if __name__ == "__main__":
    raise SystemExit(main())
