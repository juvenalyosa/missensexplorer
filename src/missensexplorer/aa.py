"""Amino-acid tables and physicochemical descriptors of a substitution.

Everything here is computed locally (no network), so it is instant and always
available, even when a remote predictor has no score for a variant.
"""

from __future__ import annotations

import math

AA1_TO_3 = {
    "A": "Ala", "R": "Arg", "N": "Asn", "D": "Asp", "C": "Cys",
    "Q": "Gln", "E": "Glu", "G": "Gly", "H": "His", "I": "Ile",
    "L": "Leu", "K": "Lys", "M": "Met", "F": "Phe", "P": "Pro",
    "S": "Ser", "T": "Thr", "W": "Trp", "Y": "Tyr", "V": "Val",
    "U": "Sec", "O": "Pyl", "*": "Ter", "X": "Xaa",
}
AA3_TO_1 = {v.upper(): k for k, v in AA1_TO_3.items()}
AA3_TO_1["STOP"] = "*"

AA_NAMES = {
    "A": "Alanine", "R": "Arginine", "N": "Asparagine", "D": "Aspartate",
    "C": "Cysteine", "Q": "Glutamine", "E": "Glutamate", "G": "Glycine",
    "H": "Histidine", "I": "Isoleucine", "L": "Leucine", "K": "Lysine",
    "M": "Methionine", "F": "Phenylalanine", "P": "Proline", "S": "Serine",
    "T": "Threonine", "W": "Tryptophan", "Y": "Tyrosine", "V": "Valine",
    "*": "Stop",
}

# Kyte & Doolittle (1982) hydropathy index.
HYDROPATHY = {
    "A": 1.8, "R": -4.5, "N": -3.5, "D": -3.5, "C": 2.5, "Q": -3.5, "E": -3.5,
    "G": -0.4, "H": -3.2, "I": 4.5, "L": 3.8, "K": -3.9, "M": 1.9, "F": 2.8,
    "P": -1.6, "S": -0.8, "T": -0.7, "W": -0.9, "Y": -1.3, "V": 4.2,
}

# Residue volumes in cubic angstroms (Zamyatnin 1972).
VOLUME = {
    "A": 88.6, "R": 173.4, "N": 114.1, "D": 111.1, "C": 108.5, "Q": 143.8,
    "E": 138.4, "G": 60.1, "H": 153.2, "I": 166.7, "L": 166.7, "K": 168.6,
    "M": 162.9, "F": 189.9, "P": 112.7, "S": 89.0, "T": 116.1, "W": 227.8,
    "Y": 193.6, "V": 140.0,
}

# Side-chain charge at physiological pH (His treated as neutral).
CHARGE = {aa: 0 for aa in HYDROPATHY}
CHARGE.update({"D": -1, "E": -1, "K": 1, "R": 1})

AA_CLASS = {
    **{aa: "hydrophobic" for aa in "AVLIM"},
    **{aa: "aromatic" for aa in "FWY"},
    **{aa: "polar" for aa in "STNQC"},
    **{aa: "positive" for aa in "KRH"},
    **{aa: "negative" for aa in "DE"},
    "G": "glycine",
    "P": "proline",
}

# Grantham (1974) composition, polarity and volume. The distance matrix is
# derived from these with Grantham's weights, which reproduces the published
# table to within rounding.
_GRANTHAM_CPV = {
    "S": (1.42, 9.2, 32.0), "R": (0.65, 10.5, 124.0), "L": (0.0, 4.9, 111.0),
    "P": (0.39, 8.0, 32.5), "T": (0.71, 8.6, 61.0), "A": (0.0, 8.1, 31.0),
    "V": (0.0, 5.9, 84.0), "G": (0.74, 9.0, 3.0), "I": (0.0, 5.2, 111.0),
    "F": (0.0, 5.2, 132.0), "Y": (0.20, 6.2, 136.0), "C": (2.75, 5.5, 55.0),
    "H": (0.58, 10.4, 96.0), "Q": (0.89, 10.5, 85.0), "N": (1.33, 11.6, 56.0),
    "K": (0.33, 11.3, 119.0), "D": (1.38, 13.0, 54.0), "E": (0.92, 12.3, 83.0),
    "M": (0.0, 5.7, 105.0), "W": (0.13, 5.4, 170.0),
}
_G_ALPHA, _G_BETA, _G_GAMMA, _G_RHO = 1.833, 0.1018, 0.000399, 50.723

_BLOSUM62_ORDER = "ARNDCQEGHILKMFPSTWYV"
_BLOSUM62_ROWS = """
 4 -1 -2 -2  0 -1 -1  0 -2 -1 -1 -1 -1 -2 -1  1  0 -3 -2  0
-1  5  0 -2 -3  1  0 -2  0 -3 -2  2 -1 -3 -2 -1 -1 -3 -2 -3
-2  0  6  1 -3  0  0  0  1 -3 -3  0 -2 -3 -2  1  0 -4 -2 -3
-2 -2  1  6 -3  0  2 -1 -1 -3 -4 -1 -3 -3 -1  0 -1 -4 -3 -3
 0 -3 -3 -3  9 -3 -4 -3 -3 -1 -1 -3 -1 -2 -3 -1 -1 -2 -2 -1
-1  1  0  0 -3  5  2 -2  0 -3 -2  1  0 -3 -1  0 -1 -2 -1 -2
-1  0  0  2 -4  2  5 -2  0 -3 -3  1 -2 -3 -1  0 -1 -3 -2 -2
 0 -2  0 -1 -3 -2 -2  6 -2 -4 -4 -2 -3 -3 -2  0 -2 -2 -3 -3
-2  0  1 -1 -3  0  0 -2  8 -3 -3 -1 -2 -1 -2 -1 -2 -2  2 -3
-1 -3 -3 -3 -1 -3 -3 -4 -3  4  2 -3  1  0 -3 -2 -1 -3 -1  3
-1 -2 -3 -4 -1 -2 -3 -4 -3  2  4 -2  2  0 -3 -2 -1 -2 -1  1
-1  2  0 -1 -3  1  1 -2 -1 -3 -2  5 -1 -3 -1  0 -1 -3 -2 -2
-1 -1 -2 -3 -1  0 -2 -3 -2  1  2 -1  5  0 -2 -1 -1 -1 -1  1
-2 -3 -3 -3 -2 -3 -3 -3 -1  0  0 -3  0  6 -4 -2 -2  1  3 -1
-1 -2 -2 -1 -3 -1 -1 -2 -2 -3 -3 -1 -2 -4  7 -1 -1 -4 -3 -2
 1 -1  1  0 -1  0  0  0 -1 -2 -2  0 -1 -2 -1  4  1 -3 -2 -2
 0 -1  0 -1 -1 -1 -1 -2 -2 -1 -1 -1 -1 -2 -1  1  5 -2 -2  0
-3 -3 -4 -4 -2 -2 -3 -2 -2 -3 -2 -3 -1  1 -4 -3 -2 11  2 -3
-2 -2 -2 -3 -2 -1 -2 -3  2 -1 -1 -2 -1  3 -3 -2 -2  2  7 -1
 0 -3 -3 -3 -1 -2 -2 -3 -3  3  1 -2  1 -1 -2 -2  0 -3 -1  4
"""
BLOSUM62: dict[tuple[str, str], int] = {}
for _a, _row in zip(_BLOSUM62_ORDER, _BLOSUM62_ROWS.strip().splitlines()):
    for _b, _v in zip(_BLOSUM62_ORDER, _row.split()):
        BLOSUM62[(_a, _b)] = int(_v)


def to_one_letter(code: str) -> str | None:
    """'Val' / 'V' / 'Ter' / '*' -> one-letter code, or None if unknown."""
    code = code.strip()
    if len(code) == 1:
        c = code.upper()
        return c if c in AA1_TO_3 else None
    return AA3_TO_1.get(code.upper())


def to_three_letter(code: str) -> str | None:
    one = to_one_letter(code)
    return AA1_TO_3.get(one) if one else None


def grantham(a: str, b: str) -> int | None:
    """Grantham distance between two residues (0 = identical, 215 = max)."""
    if a not in _GRANTHAM_CPV or b not in _GRANTHAM_CPV:
        return None
    ca, pa, va = _GRANTHAM_CPV[a]
    cb, pb, vb = _GRANTHAM_CPV[b]
    d = math.sqrt(
        _G_ALPHA * (ca - cb) ** 2 + _G_BETA * (pa - pb) ** 2 + _G_GAMMA * (va - vb) ** 2
    )
    return round(_G_RHO * d)


def grantham_class(distance: int | None) -> str | None:
    """Li et al. (1984) categories."""
    if distance is None:
        return None
    if distance <= 50:
        return "conservative"
    if distance <= 100:
        return "moderately conservative"
    if distance <= 150:
        return "moderately radical"
    return "radical"


def substitution_profile(ref: str, alt: str) -> dict:
    """Physicochemical description of a single amino-acid substitution."""
    out: dict = {
        "grantham": None, "grantham_class": None, "blosum62": None,
        "hydropathy_change": None, "volume_change": None, "charge_change": None,
        "aa_class_change": None, "physchem_notes": None,
    }
    if ref not in HYDROPATHY or alt not in HYDROPATHY or ref == alt:
        return out
    g = grantham(ref, alt)
    out["grantham"] = g
    out["grantham_class"] = grantham_class(g)
    out["blosum62"] = BLOSUM62.get((ref, alt))
    out["hydropathy_change"] = round(HYDROPATHY[alt] - HYDROPATHY[ref], 2)
    out["volume_change"] = round(VOLUME[alt] - VOLUME[ref], 1)
    out["charge_change"] = CHARGE[alt] - CHARGE[ref]
    out["aa_class_change"] = f"{AA_CLASS[ref]} -> {AA_CLASS[alt]}"

    notes = []
    if CHARGE[ref] != CHARGE[alt]:
        notes.append(f"charge {CHARGE[ref]:+d} -> {CHARGE[alt]:+d}")
    if abs(out["hydropathy_change"]) >= 3:
        notes.append("large hydropathy shift" + (" (more hydrophobic)" if out["hydropathy_change"] > 0 else " (more polar)"))
    if abs(out["volume_change"]) >= 50:
        notes.append("large size change" + (" (bulkier)" if out["volume_change"] > 0 else " (smaller)"))
    if alt == "P":
        notes.append("proline introduced (breaks helices/strands)")
    if ref == "P":
        notes.append("proline lost")
    if ref == "G":
        notes.append("glycine lost (backbone flexibility)")
    if alt == "G":
        notes.append("glycine introduced")
    if ref == "C":
        notes.append("cysteine lost (possible disulfide/metal site)")
    if alt == "C":
        notes.append("cysteine introduced (unpaired thiol)")
    out["physchem_notes"] = "; ".join(notes) or None
    return out
