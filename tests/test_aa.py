import pytest

from missensexplorer.aa import BLOSUM62, grantham, grantham_class, substitution_profile, to_one_letter


@pytest.mark.parametrize("a,b,expected", [
    ("C", "W", 215), ("L", "I", 5), ("S", "R", 110), ("D", "E", 45), ("F", "Y", 22),
    ("G", "W", 184), ("V", "E", 121), ("R", "C", 180), ("R", "H", 29),
])
def test_grantham_matches_published_table(a, b, expected):
    assert abs(grantham(a, b) - expected) <= 1
    assert grantham(a, b) == grantham(b, a)


def test_blosum62_symmetric_and_known_values():
    assert len(BLOSUM62) == 400
    assert all(BLOSUM62[(a, b)] == BLOSUM62[(b, a)] for a, b in BLOSUM62)
    assert BLOSUM62[("W", "W")] == 11 and BLOSUM62[("V", "E")] == -2 and BLOSUM62[("I", "V")] == 3


def test_substitution_profile_flags():
    p = substitution_profile("V", "E")
    assert p["charge_change"] == -1 and p["grantham_class"] == "moderately radical"
    assert "charge" in p["physchem_notes"]
    assert "proline introduced" in substitution_profile("L", "P")["physchem_notes"]
    assert "glycine lost" in substitution_profile("G", "D")["physchem_notes"]
    assert substitution_profile("V", "V")["grantham"] is None


def test_codes():
    assert to_one_letter("Val") == "V" and to_one_letter("ter") == "*" and to_one_letter("q") == "Q"
    assert grantham_class(10) == "conservative" and grantham_class(200) == "radical"
