"""Offscreen smoke tests for the Qt desktop app (skipped when PySide6 is not installed)."""

import os
import tempfile

import pandas as pd
import pytest

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
os.environ.setdefault("MISSENSEXPLORER_NO_WEBENGINE", "1")
QtWidgets = pytest.importorskip("PySide6.QtWidgets")

# Keep the tests away from the user's real preferences: the window reads them and saves them on close.
os.environ["MISSENSEXPLORER_SETTINGS"] = os.path.join(tempfile.mkdtemp(prefix="mx-settings-"), "settings.ini")

from missensexplorer.gui import i18n  # noqa: E402
from missensexplorer.predict import COLUMNS, PredictionResult  # noqa: E402


@pytest.fixture(scope="module")
def app():
    return QtWidgets.QApplication.instance() or QtWidgets.QApplication([])


def _result():
    rows = [
        {"input": "BRAF V600E", "gene": "BRAF", "protein_change": "V600E", "hgvsp": "p.Val600Glu",
         "consequence": "missense_variant", "prediction": "Likely damaging", "effect_score": 1.0,
         "prediction_basis": "4 damaging / 0 benign predictors (weighted 100% damaging) | AlphaMissense 0.99",
         "alphamissense": 0.99, "alphamissense_class": "likely_pathogenic", "revel": 0.93, "cadd_phred": 29.8,
         "ref_aa": "V", "alt_aa": "E", "grantham": 121, "grantham_class": "moderately radical", "blosum62": -2,
         "hydropathy_change": -7.7, "volume_change": 1.5, "charge_change": -1, "aa_class_change": "hydrophobic -> negative",
         "protein_position": 600, "protein_length": 766, "protein_window": "GDFGLAT[V>E]KSRWSGS", "uniprot": "P15056",
         "clinvar": "Pathogenic", "clinvar_stars": 2, "gnomad_af": 4e-6, "chrom": "7", "pos": 140753336,
         "ref": "A", "alt": "T", "acmg_pp3_bp4": "PP3_Moderate (REVEL 0.930)", "acmg_points": 3,
         "acmg_ps1_pm5": "PM5 (p.Val600Gly: Pathogenic, 3★, VCV000040389)",
         "acmg_codes": "PP3_Moderate +3, PM5 +2", "population_frequency": "rare (grpmax 0.0004%)"},
        {"input": "TP53 P72R", "gene": "TP53", "protein_change": "P72R", "consequence": "missense_variant",
         "prediction": "Likely benign", "effect_score": 0.0, "protein_position": 72, "protein_length": 393},
        {"input": "TP53 R213*", "gene": "TP53", "protein_change": "R213*", "consequence": "stop_gained",
         "prediction": "Loss of function (truncating)", "effect_score": 1.0, "protein_position": 213,
         "protein_length": 393, "ref_aa": "R", "alt_aa": "*"},
    ]
    table = pd.DataFrame(rows, columns=COLUMNS)
    errors = pd.DataFrame([{"input": "FOO X1Y", "source": "text", "error": "unknown gene or transcript symbol"}])
    wt = "M" * 212 + "R" + "A" * 180
    return PredictionResult(table=table, errors=errors, assembly="GRCh38", stats={"inputs": 4, "rows": 3,
                            "likely_damaging": 2, "errors": 1, "seconds": 1.2},
                            sequences={2: {"wt": wt, "mut": wt[:212], "note": "truncated", "header": "TP53|R213*"}})


@pytest.mark.parametrize("lang", ["es", "en"])
def test_window_shows_results_and_details(app, lang):
    i18n.set_language(lang)
    from missensexplorer.gui.window import MainWindow

    w = MainWindow()
    w.show_result(_result())
    assert w.stack.currentIndex() == 1
    assert w.proxy.rowCount() == 3
    assert w.k_cv.value.text() == "1"
    assert set(w.segments) == {"Damaging", "Benign / silent"}

    w.segments["Damaging"].setChecked(False)          # filter by status
    assert w.proxy.rowCount() == 1
    w.segments["Damaging"].setChecked(True)
    w.search.setText("tp53")                           # free-text search
    assert w.proxy.rowCount() == 2
    w.search.clear()

    w.select_index(2)                                  # truncating variant: sequence tab filled
    assert w.current_idx == 2
    assert "212 aa" in w.seq_view.toPlainText()
    w.select_index(0)
    texts = " ".join(lbl.text() for lbl in w.detail.findChildren(QtWidgets.QLabel))
    assert "BRAF p.Val600Glu" in texts
    assert ("Probablemente dañina" if lang == "es" else "Likely damaging") in texts
    assert "PP3_Moderate" in texts and ("rara (grpmax" if lang == "es" else "rare (grpmax") in texts
    assert ("PM5 (p.Val600Gly: patogénica, 3★" if lang == "es" else "PM5 (p.Val600Gly: Pathogenic, 3★") in texts
    assert "PP3_Moderate +3, PM5 +2" in texts
    from missensexplorer.gui.widgets import ContributionChart
    charts = w.detail.findChildren(ContributionChart)
    assert len(charts) == 1 and {n for n, _ in charts[0].items} == {"AlphaMissense", "REVEL", "CADD"}
    assert len(w.lolli.points) == 1 and w.lolli.selected == 0
    assert not w.grab().isNull()
    w.close()


def test_language_switch_keeps_state(app):
    from missensexplorer.gui.window import MainWindow

    i18n.set_language("es")
    old = MainWindow()
    old.text.setPlainText("BRAF V600E")
    old.show_result(_result())
    old.select_index(1)
    i18n.set_language("en")
    new = MainWindow()
    new.transfer_state(old)
    assert new.text.toPlainText() == "BRAF V600E"
    assert new.current_idx == 1 and new.proxy.rowCount() == 3
    assert new.title.text() == "Results"
    old.close()
    new.close()
