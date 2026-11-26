"""ClinGen expert panel (VCEP) specifications of the ACMG/AMP criteria, by gene.

ClinGen's Variant Curation Expert Panels specify the ACMG/AMP criteria for their genes, and Pejaver et
al. (2022) note that laboratories could follow a panel's recommended tool and thresholds where one exists.
The panels' thresholds are free text; what is used here comes from ClinGen's Criteria Specification
Registry (CSpec; vcep_params.py, built by benchmarks/make_vcep_table.py):

- the strengths each panel allows for each criterion, and the PP3/BP4 it does not use for missense
  variants although CSpec marks them applicable (the text says so);
- its BA1 allele-frequency threshold, read from the text;
- its conditions the tool cannot check, quoted;
- for the ten panels whose PP3/BP4, PM2 and BA1 rules were encoded by hand (vcep_rules.py), what those
  rules give for the variant.

For a variant in a gene with a released specification, the output names the specification (with its
disease and link) and gives one note per specification with every point where the tool's evidence differs
from what the panel allows. The tool's PP3/BP4 is its own calibration, not the panel's rule: the note says
so whenever the tool gives PP3 or BP4.
"""

from __future__ import annotations

import re
from numbers import Real
from typing import Any, cast

from . import vcep_rules
from ._types import Row
from .vcep_params import VCEP as _VCEP

VCEP = cast(dict[str, Any], _VCEP)

POINTS = {"Supporting": 1, "Moderate": 2, "Strong": 4, "Very Strong": 8, "Stand Alone": 8}
NAMES = {1: "Supporting", 2: "Moderate", 3: "Moderate+", 4: "Strong", 8: "Very Strong"}
CODE_POINTS = {
    "PS1": 4,
    "PS1_Moderate": 2,
    "PM5": 2,
    "PM5_Supporting": 1,
    "PM2_Supporting": 1,
}
SEP = " | "  # between the notes of different specifications in the acmg_notes column
MAX_QUOTES = 4  # conditions quoted per criterion; the rest are counted
# BA1 texts with conditions beyond the threshold: allele counts, sex, other datasets, exclusions, founders
BA1_CONDITIONS = re.compile(
    r"alleles|allele count|\bmales?\b|hemizyg|homozyg|Exome Sequencing Project|\bESP\b|ExAC|1000 Genomes|"
    r"exclu|founder|except|at least|NOT MET|fail",
    re.I,
)


def specifications(gene: str | None) -> list[dict[str, Any]]:
    """The released specifications that cover `gene`, once each (empty when none does)."""
    seen: set[str] = set()
    out = []
    for s in VCEP["genes"].get(str(gene or ""), []):
        if s["id"] not in seen:
            seen.add(s["id"])
            out.append(s)
    return out


def _short(panel: str) -> str:
    return " ".join(panel.replace("Variant Curation Expert Panel", "VCEP").split())


def _diseases(s: dict[str, Any], n: int) -> str:
    """The specification's diseases, with 'X type 2A', 'X type 2B' as 'X types 2A, 2B', cut to `n` characters."""
    order: list[str] = []
    types: dict[str, list[str]] = {}
    for d in s.get("diseases") or []:
        m = re.match(r"^(.*\btype) (\S+)$", d)
        key = m[1] if m else d
        if key not in types:
            order.append(key)
            types[key] = []
        if m:
            types[key].append(m[2])
    parts = [
        f"{k}s {', '.join(types[k])}" if len(types[k]) > 1 else f"{k} {types[k][0]}" if types[k] else k for k in order
    ]
    return _cut(", ".join(parts), n)


def _name(s: dict[str, Any], several: bool) -> str:
    """'Congenital Myopathies VCEP v2.0' plus, when the gene has several specifications, what tells them apart."""
    name = f"{_short(s['panel'])} v{s['version']}"
    if several:
        what = ", ".join(x for x in (_diseases(s, 60), s.get("inheritance")) if x)
        name += f" ({what}; {s['id']})" if what else f" ({s['id']})"
    return name


def _cut(text: str, n: int) -> str:
    """`text` cut to at most `n` characters at a word boundary."""
    return text if len(text) <= n else text[:n].rsplit(" ", 1)[0].rstrip(",;") + "…"


def describe(gene: str | None) -> str | None:
    """'Familial Hypercholesterolemia VCEP v1.2, hypercholesterolemia, familial, 1 (https://...)' per specification."""
    specs = specifications(gene)
    out = []
    for s in specs:
        dis = _diseases(s, 80)
        inh = f", {s['inheritance']}" if s.get("inheritance") else ""
        out.append(f"{_short(s['panel'])} v{s['version']}{', ' + dis if dis else ''}{inh} ({s['url']})")
    return "; ".join(out) or None


def _tool_evidence(row: Row) -> list[tuple[str, str, int]]:
    """(criterion, code as written, points) of the ACMG criteria the tool gave this row."""
    out = []
    pts = _num(row.get("acmg_points"))
    if pts is not None and int(pts):
        p = int(pts)
        crit = "PP3" if p > 0 else "BP4"
        out.append((crit, f"{crit}_{NAMES.get(abs(p), str(abs(p)))}", abs(p)))
    for col in ("acmg_ps1_pm5", "acmg_pm2"):
        code = str(row.get(col) or "").split(" ")[0]
        if code in CODE_POINTS:
            out.append((code.split("_")[0], code, CODE_POINTS[code]))
    return out


def _num(x: Any) -> float | None:
    """A finite number (NumPy numbers too), or None."""
    if isinstance(x, bool) or not isinstance(x, Real):
        return None
    v = float(x)
    return v if v == v else None


def _pct(x: float) -> str:
    return f"{100 * x:.3g}%"


def _ba1(s: dict[str, Any], row: Row) -> str | None:
    """Where the panel's BA1 threshold and the tool's BA1 (grpmax above 5%) disagree."""
    from .acmg import ba1_grpmax  # the grpmax the tool's own BA1 uses (groups with >= 2,000 alleles when known)

    b = s.get("ba1")
    g = _num(ba1_grpmax(row))
    if g is None:
        g = _num(row.get("gnomad_popmax_af"))
    if not b or g is None:
        return None
    tool = str(row.get("population_frequency") or "").startswith("common: BA1")
    if b.get("value") is None:
        return f"its own BA1 rule: “{_cut(b['quote'], 120)}”" if tool or g >= 0.001 else None
    t, op = float(b["value"]), b.get("op") or ">"
    met = g >= t if op == ">=" else g > t
    if tool and not met:
        return f"BA1 only {op} {_pct(t)}: grpmax {_pct(g)} does not meet it"
    if met and not tool and b.get("faf"):
        # the grpmax allele frequency is at least the filtering allele frequency: "met" is not established
        return (
            f"grpmax allele frequency {_pct(g)} is above its BA1 threshold ({op} {_pct(t)}), which it applies to the "
            "filtering allele frequency (not computed here): check the filtering allele frequency before BA1"
        )
    if met and not tool:
        text = f"grpmax {_pct(g)} is at or above its BA1 threshold ({op} {_pct(t)}): check its other BA1 conditions"
        if BA1_CONDITIONS.search(b.get("quote") or ""):
            text += f" (“{_cut(b['quote'], 200)}”)"
        return text
    return None


def ba1_above_5pct(gene: str | None) -> float | None:
    """The BA1 threshold when every specification of the gene sets it above 5% (VWF, ABCA4), else None.

    BA1 at 5% applies to genes without a gene-specific BA1 rule (Ghosh et al. 2018); these panels set it
    higher, so a variant between 5% and their threshold does not meet BA1.
    """
    vals = [_num((s.get("ba1") or {}).get("value")) for s in specifications(gene)]
    nums = [v for v in vals if v is not None]
    if vals and len(nums) == len(vals) and all(v > 0.05 for v in nums):
        return min(nums)
    return None


def note(s: dict[str, Any], row: Row, several: bool = False) -> str | None:
    """One specification's note on the row: where the tool's evidence differs from what the panel allows."""
    parts: list[str] = []
    tool = _tool_evidence(row)
    for crit, code, pts in tool:
        if crit in s.get("not_missense", {}):
            parts.append(f"does not use {crit} for missense variants (“{_cut(s['not_missense'][crit], 100)}”)")
            continue
        allowed = s["allowed"].get(crit)
        if allowed is None:
            continue
        if not allowed:
            parts.append(f"does not use {crit}")
            continue
        levels = sorted(POINTS[a] for a in allowed)
        if pts > levels[-1]:
            parts.append(f"allows {crit} up to {NAMES[levels[-1]]}: the tool's {code} ({pts} points) is stronger")
        elif pts < levels[0] and crit == "PM2":
            # confirmed absence meets a moderate PM2 threshold
            parts.append(f"uses PM2 only at {', '.join(allowed)}: the tool's {code} is weaker")
        elif pts < levels[0] and crit in ("PS1", "PM5"):
            # the tool lowered the strength because the reference is likely pathogenic or has one review star
            parts.append(
                f"has no {crit} below {allowed[0]}: the tool lowered {crit} because the ClinVar reference is likely "
                "pathogenic or has one review star, and the panel judges the reference under its own specification"
            )
        elif pts < levels[0]:
            parts.append(f"uses {crit} only at {', '.join(allowed)} (its own thresholds)")
        elif pts not in levels:
            parts.append(f"has no {pts}-point {crit} ({', '.join(allowed)}): the tool's {code} is not one of them")
    rule = vcep_rules.RULES.get(str(row.get("gene") or ""))
    if rule and rule.get("spec") in (s["id"], "RASopathy") and str(s["panel"]).startswith(rule["panel"]):
        if got := vcep_rules.describe(row):
            parts.append(f"its rules give {got}")
    elif any(
        c in ("PP3", "BP4") and c not in s.get("not_missense", {}) and s["allowed"].get(c) != [] for c, _, _ in tool
    ):
        parts.append("it sets its own predictor and thresholds for PP3/BP4 (the tool's PP3/BP4 is its own calibration)")
    if b := _ba1(s, row):
        parts.append(b)
    for crit, _code, _ in tool:
        quotes = s.get("conditions", {}).get(crit, [])
        for q in quotes[:MAX_QUOTES]:
            parts.append(f"on {crit}: “{q}”")
        if len(quotes) > MAX_QUOTES:
            parts.append(f"on {crit}: {len(quotes) - MAX_QUOTES} more conditions in the specification")
    if not parts:
        return None
    return f"{_name(s, several)}: {'; '.join(parts)} ({s['url']})"


def notes(row: Row) -> list[str]:
    """One note per specification covering the row's gene, on what differs from the tool's evidence."""
    specs = specifications(row.get("gene"))
    return [n for s in specs if (n := note(s, row, len(specs) > 1))]
