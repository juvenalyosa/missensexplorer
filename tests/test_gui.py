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
        {
            "input": "BRAF V600E",
            "gene": "BRAF",
            "protein_change": "V600E",
            "hgvsp": "p.Val600Glu",
            "consequence": "missense_variant",
            "prediction": "Likely damaging",
            "effect_score": 1.0,
            "prediction_basis": "4 damaging / 0 benign predictors (weighted 100% damaging) | AlphaMissense 0.99",
            "alphamissense": 0.99,
            "alphamissense_class": "likely_pathogenic",
            "revel": 0.93,
            "cadd_phred": 29.8,
            "ref_aa": "V",
            "alt_aa": "E",
            "grantham": 121,
            "grantham_class": "moderately radical",
            "blosum62": -2,
            "hydropathy_change": -7.7,
            "volume_change": 1.5,
            "charge_change": -1,
            "aa_class_change": "hydrophobic -> negative",
            "protein_position": 600,
            "protein_length": 766,
            "protein_window": "GDFGLAT[V>E]KSRWSGS",
            "uniprot": "P15056",
            "clinvar": "Pathogenic",
            "clinvar_stars": 2,
            "gnomad_af": 4e-6,
            "chrom": "7",
            "pos": 140753336,
            "ref": "A",
            "alt": "T",
            "acmg_pp3_bp4": "PP3_Moderate (REVEL 0.930)",
            "acmg_points": 3,
            "acmg_ps1_pm5": "PM5 (p.Val600Gly: Pathogenic, 3★, VCV000040389)",
            "acmg_codes": "PP3_Moderate +3, PM5 +2",
            "population_frequency": "rare (grpmax 0.0004%)",
            "score_with_frequency": 0.97,
        },
        {
            "input": "TP53 P72R",
            "gene": "TP53",
            "protein_change": "P72R",
            "consequence": "missense_variant",
            "prediction": "Likely benign",
            "effect_score": 0.0,
            "protein_position": 72,
            "protein_length": 393,
        },
        {
            "input": "TP53 R213*",
            "gene": "TP53",
            "protein_change": "R213*",
            "consequence": "stop_gained",
            "prediction": "Loss of function (truncating)",
            "effect_score": 1.0,
            "protein_position": 213,
            "protein_length": 393,
            "ref_aa": "R",
            "alt_aa": "*",
        },
    ]
    table = pd.DataFrame(rows, columns=COLUMNS)
    errors = pd.DataFrame([{"input": "FOO X1Y", "source": "text", "error": "unknown gene or transcript symbol"}])
    wt = "M" * 212 + "R" + "A" * 180
    return PredictionResult(
        table=table,
        errors=errors,
        assembly="GRCh38",
        stats={"inputs": 4, "rows": 3, "likely_damaging": 2, "errors": 1, "seconds": 1.2},
        sequences={2: {"wt": wt, "mut": wt[:212], "note": "truncated", "header": "TP53|R213*"}},
    )


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

    w.segments["Damaging"].setChecked(False)  # filter by status
    assert w.proxy.rowCount() == 1
    w.segments["Damaging"].setChecked(True)
    w.search.setText("tp53")  # free-text search
    assert w.proxy.rowCount() == 2
    w.search.clear()

    w.select_index(2)  # truncating variant: sequence tab filled
    assert w.current_idx == 2
    assert "212 aa" in w.seq_view.toPlainText()
    w.select_index(0)
    texts = " ".join(lbl.text() for lbl in w.detail.findChildren(QtWidgets.QLabel))
    assert "BRAF p.Val600Glu" in texts
    assert ("Probablemente dañina" if lang == "es" else "Likely damaging") in texts
    assert "PP3_Moderate" in texts
    assert ("rara (grpmax" if lang == "es" else "rare (grpmax") in texts
    assert ("PM5 (p.Val600Gly: patogénica, 3★" if lang == "es" else "PM5 (p.Val600Gly: Pathogenic, 3★") in texts
    assert "PP3_Moderate +3, PM5 +2" in texts
    assert ("con frecuencia 0.97" if lang == "es" else "with frequency 0.97") in texts
    from missensexplorer.gui.widgets import ContributionChart

    charts = w.detail.findChildren(ContributionChart)
    assert len(charts) == 1
    assert {n for n, _ in charts[0].items} == {"AlphaMissense", "REVEL", "CADD"}
    assert len(w.lolli.points) == 1
    assert w.lolli.selected == 0
    assert not w.grab().isNull()
    w.close()


def test_acmg_note_names_pp3_bp4_thresholds_only_when_given(app):
    i18n.set_language("en")
    from missensexplorer.gui.window import MainWindow

    res = _result()
    res.table.loc[2, "population_frequency"] = "rare (grpmax 0.00009%)"  # stop-gained: ACMG evidence, no PP3/BP4
    w = MainWindow()
    w.show_result(res)

    def detail_text() -> str:
        return " ".join(lbl.text() for lbl in w.detail.findChildren(QtWidgets.QLabel))

    w.select_index(0)  # PP3_Moderate (REVEL 0.930)
    assert i18n.t("acmg_note") in detail_text()
    w.select_index(2)
    text = detail_text()
    assert "rare (grpmax" in text
    assert i18n.t("acmg_note_general") in text
    assert "REVEL" not in text
    w.close()


@pytest.mark.parametrize("lang", ["es", "en"])
def test_skipped_note_counts_records_that_failed_filter(app, lang):
    i18n.set_language(lang)
    from missensexplorer.gui.window import MainWindow

    res = _result()
    res.stats.update({"file_skipped_filter": 1, "file_skipped_not_carried": 2})
    w = MainWindow()
    w.show_result(res)
    note = w.skipped.text()
    assert (
        "2 no portadas por la muestra · 1 sin FILTER = PASS"
        if lang == "es"
        else "2 not carried by the sample(s) · 1 not passing FILTER"
    ) in note
    w.close()


@pytest.mark.parametrize(
    "text",
    [
        "the amino acid at this position is not the one given (check the position and the isoform)",
        "position outside the protein (check the position and the isoform)",
        "a frameshift written as a protein change cannot be traced to one DNA change: give the c. change",
        "2 nucleotide changes give this protein change; each is reported",
        "protein notation: other DNA changes can give the same F508del; this is the one Ensembl maps "
        "(c.1522_1524del). Give the c. or genomic change to analyse a specific one",
        "not checked (no MANE Select transcript, on which ClinVar names its variants)",
        "protein change given as p.Gly396Asp (another isoform's numbering); on the reported transcript it is G368D",
        "protein change given as p.Ala175His exists only on another isoform; on the reported transcript this DNA "
        "change is C182S: check the input",
    ],
)
def test_new_pipeline_messages_are_translated(text):
    i18n.set_language("es")
    out = i18n.tr_text(text)
    assert out != text
    assert not any(w in out for w in ("the ", "position", "change", "transcript,")), out
    i18n.set_language("en")


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
    assert new.current_idx == 1
    assert new.proxy.rowCount() == 3
    assert new.title.text() == "Results"
    old.close()
    new.close()


def test_prediction_column_sorts_rows_without_a_score(app):
    from missensexplorer.gui.models import COLUMNS as GUI_COLUMNS
    from missensexplorer.gui.models import SORT_ROLE, ResultsModel

    df = pd.DataFrame(
        [
            {"prediction": "Likely damaging", "effect_score": 0.9},
            {"prediction": "Likely damaging", "effect_score": float("nan")},  # missing score
        ],
        columns=COLUMNS,
    )
    model = ResultsModel()
    model.set_frame(df)
    col = [c[0] for c in GUI_COLUMNS].index("prediction")
    keys = [model.data(model.index(i, col), SORT_ROLE) for i in range(2)]
    assert all(k == k for k in keys)  # no NaN: a NaN key makes sorting inconsistent
    assert keys[0] < keys[1]  # the scored variant sorts first within its status
