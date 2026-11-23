"""ACMG PS1 and PM5 from ClinVar: pathogenic missense variants at the same amino-acid residue.

PS1  the same amino-acid change as an established pathogenic variant, made by a different nucleotide
     change (Strong, 4 points; PS1_Moderate, 2 points, when that variant is only likely pathogenic);
PM5  a different missense change at a residue where a missense change is established as pathogenic
     (Moderate, 2 points; PM5_Supporting, 1 point, when only likely pathogenic ones are known).
The lower strengths for likely pathogenic reference variants follow several ClinGen expert panel
specifications (e.g. glaucoma, HNF4A). Only one of the two is given (the stronger), never both.

Reference variants are ClinVar germline classifications Pathogenic, Pathogenic/Likely pathogenic or
Likely pathogenic with review criteria (>= 1 star, no conflicts). They are matched on the reported
transcript (MANE Select RefSeq, as in ClinVar's variant names), the protein position and the
reference amino acid; the variant itself never counts. When the variant may act on splicing
(SpliceAI >= 0.2), PS1/PM5 are not given: its effect need not be the amino-acid change.

Not checked automatically: whether the reference variant's classification is independent of PS1/PM5,
and whether the reference variant itself acts on splicing (check it for codons at exon edges).
"""

from __future__ import annotations

import re
from dataclasses import dataclass

from .aa import to_one_letter

TITLE = re.compile(r"^(N[MR]_\d+)(?:\.\d+)?\([^)]*\):(c\.\S+) \(p\.([A-Z][a-z]{2})(\d+)([A-Z][a-z]{2})\)$")
PATHOGENIC = {"Pathogenic": "P", "Pathogenic/Likely pathogenic": "P", "Likely pathogenic": "LP"}
REVIEWED = {"practice guideline": 4, "reviewed by expert panel": 3,
            "criteria provided, multiple submitters, no conflicts": 2, "criteria provided, single submitter": 1}
SPLICE_CAUTION = 0.2
POINTS = {"PS1": 4, "PS1_Moderate": 2, "PM5": 2, "PM5_Supporting": 1}
WINDOW = 2          # bases on each side of the variant: covers its whole codon on either strand


@dataclass
class Reference:
    accession: str
    cdna: str             # c. change of the reference variant
    protein: str          # p.Arg156His
    level: str            # P or LP
    stars: int
    same_change: bool     # True: same amino-acid change (PS1); False: another change here (PM5)


def window(row: dict) -> tuple[str, int, int] | None:
    """Genomic window to search for a missense SNV, or None when PS1/PM5 do not apply to the row."""
    if "missense_variant" not in str(row.get("consequence") or "") or not row.get("protein_position"):
        return None
    ref, alt, pos = str(row.get("ref") or ""), str(row.get("alt") or ""), row.get("pos")
    if len(ref) != 1 or len(alt) != 1 or not pos:
        return None
    return str(row.get("chrom")), int(pos) - WINDOW, int(pos) + WINDOW


def references(row: dict, records: dict[str, dict], assembly: str) -> list[Reference]:
    """Pathogenic or likely pathogenic ClinVar missense variants at the row's residue."""
    nm = str(row.get("mane_select") or "").split(".")[0]
    pos, ref_aa, alt_aa = row.get("protein_position"), row.get("ref_aa"), row.get("alt_aa")
    own = (str(row.get("chrom")), int(row.get("pos") or 0), str(row.get("ref")), str(row.get("alt")))
    own_cdna = str(row.get("hgvsc") or "").split(":")[-1]   # MANE: Ensembl and RefSeq share c. numbering
    out = []
    for acc, rec in records.items():
        level = PATHOGENIC.get(rec.get("classification") or "")
        stars = REVIEWED.get(rec.get("review_status") or "", 0)
        m = TITLE.match(rec.get("title") or "")
        if not level or not stars or not m or m.group(1) != nm or int(m.group(4)) != pos:
            continue
        r_aa, a_aa = to_one_letter(m.group(3)), to_one_letter(m.group(5))
        if r_aa != ref_aa or not a_aa or a_aa in (r_aa, "*"):
            continue                                   # another residue, synonymous or nonsense
        if m.group(2) == own_cdna or _same_allele(rec, own, assembly):
            continue                                   # the variant itself
        out.append(Reference(acc, m.group(2), f"p.{m.group(3)}{m.group(4)}{m.group(5)}", level, stars,
                             a_aa == alt_aa))
    return sorted(out, key=lambda x: (not x.same_change, x.level != "P", -x.stars, x.accession))


def _same_allele(rec: dict, own: tuple, assembly: str) -> bool:
    """Whether a ClinVar record is the row's own allele, by its canonical SPDI (GRCh38 only)."""
    spdi = rec.get("spdi")
    if not spdi or assembly != "GRCh38":
        return False
    from .clinvar import _parse_spdi

    p = _parse_spdi(spdi)
    return bool(p) and (p[0], p[1] + 1, p[2], p[3]) == own


def evidence(row: dict, refs: list[Reference] | None) -> tuple[str | None, int | None]:
    """('PM5 (p.Arg156His: Pathogenic, 2★, VCV000012345)', 2), ('none (...)', 0), or (None, None)
    when PS1/PM5 do not apply to the row or were not checked (`refs` None)."""
    if refs is None:
        return None, None
    if (row.get("spliceai_max") or 0) >= SPLICE_CAUTION:
        return f"not applied (SpliceAI {row['spliceai_max']:.2f}: the variant may act on splicing)", 0
    if not refs:
        return "none (no pathogenic missense variant at this residue in ClinVar)", 0
    same = [x for x in refs if x.same_change]
    other = [x for x in refs if not x.same_change]
    options = []
    if same:
        options.append(("PS1" if same[0].level == "P" else "PS1_Moderate", same))
    if other:
        options.append(("PM5" if other[0].level == "P" else "PM5_Supporting", other))
    code, used = max(options, key=lambda o: POINTS[o[0]])     # PS1 first on a tie (the more specific)
    shown = "; ".join(_describe(x) for x in used[:2]) + (f"; and {len(used) - 2} more" if len(used) > 2 else "")
    return f"{code} ({shown})", POINTS[code]


def _describe(x: Reference) -> str:
    level = "Pathogenic" if x.level == "P" else "Likely pathogenic"
    what = f"same change by {x.cdna}" if x.same_change else x.protein
    return f"{what}: {level}, {x.stars}★, {x.accession}"
