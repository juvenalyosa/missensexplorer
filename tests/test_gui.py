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
            "acmg_pp3_bp4": "PP3_Moderate+ (REVEL 0.930)",
            "acmg_points": 3,
            "acmg_ps1_pm5": "PM5 (p.Val600Gly: Pathogenic, 3★, VCV000040389)",
            "acmg_codes": "PP3_Moderate+ +3, PM5 +2",
            "population_frequency": "rare (grpmax 0.0004%)",
            "score_with_frequency": 0.97,
            "oncogenicity": "Likely oncogenic",
            "oncogenicity_points": 6,
            "oncogenicity_codes": "OS3 +4, OP1 +1, OP4 +1",
            "oncogenicity_evidence": "OS3 (hotspot V600: 897 tumours, 833 with E; cancerhotspots.org); OP1 "
            "(calibrated computational evidence reaches PP3); OP4 (absent from gnomAD); BRAF: activating (UniProt: "
            "proto-oncogene); in an activating gene, a damage prediction does not tell activation from "
            "inactivation; not assessed automatically: functional studies beyond CIViC, OM1 functional domain, OP2, "
            "tumour type",
            "clinvar_oncogenicity": "Oncogenic",
            "clinical_evidence_level": "A",
            "clinical_evidence": "88 CIViC items (level A 6, B 55, C 16, D 11); best: level A predictive, "
            "Sensitivity/Response: Vemurafenib in Skin Melanoma (EID1409)",
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
    assert "PP3_Moderate+" in texts
    assert ("rara (grpmax" if lang == "es" else "rare (grpmax") in texts
    assert ("PM5 (p.Val600Gly: patogénica, 3★" if lang == "es" else "PM5 (p.Val600Gly: Pathogenic, 3★") in texts
    assert "PP3_Moderate+ +3, PM5 +2" in texts
    assert ("con frecuencia 0.97" if lang == "es" else "with frequency 0.97") in texts
    assert ("Probablemente oncogénica · +6 puntos" if lang == "es" else "Likely oncogenic · +6 points") in texts
    assert ("hotspot V600: 897 tumores, 833 con E" if lang == "es" else "hotspot V600: 897 tumours") in texts
    assert ("activador (UniProt: protooncogén)" if lang == "es" else "activating (UniProt: proto-oncogene)") in texts
    assert (
        "88 evidencias en CIViC (nivel A 6, B 55, C 16, D 11); la mejor: nivel A, predictiva, sensibilidad/respuesta: "
        "Vemurafenib en Skin Melanoma (EID1409)"
        if lang == "es"
        else "best: level A predictive, Sensitivity/Response: Vemurafenib in Skin Melanoma"
    ) in texts
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
        "an insertion written as a protein change cannot be traced to one DNA change: give the c. change",
        "several DNA changes give this protein change and Ensembl did not resolve them: give the c. change",
        "Ensembl places this protein change outside the coding sequence of the gene's main transcript "
        "(check the position and the isoform)",
        "2 nucleotide changes give this protein change; each is reported",
        "protein notation: other DNA changes can give the same F508del; this is the one Ensembl maps "
        "(c.1522_1524del). Give the c. or genomic change to analyse a specific one",
        "not checked (no MANE Select transcript, on which ClinVar names its variants)",
        "in an activating gene, a damage prediction does not tell activation from inactivation",
        "OP1 (calibrated computational evidence reaches PP3)",
        "OVS1_Strong (NMD escape likely (last exon); the truncation removes 24% of the protein (ClinGen LoF "
        "guidance: > 10% strong))",
        "OVS1 (canonical splice site; skipping exon 4 would shift the reading frame)",
        "OVS1_Moderate (start codon lost: another in-frame start codon may be used (ClinGen LoF guidance))",
        "OVS1 not applied (not a loss-of-function gene: not a known cancer gene)",
        "CIViC oncogenic evidence from clinical observations (not functional studies) not used",
        "OS1/OM4 not checked (ClinVar oncogenicity search unavailable)",
        "none (model 0.050; BP4_Moderate+ withheld: SpliceAI 0.30 does not rule out a splicing effect)",
        "frequent: BS1 candidate (grpmax 6.1%); BA1 not met: 4.2% at most in the gnomAD groups with >= 2,000 alleles",
        "PP3 (4 points) and PM5 may partly reflect the same residue-level information (Pejaver et al. 2022): weigh "
        "their combination with care; if PM1 is added, PP3 + PM1 should not exceed 4 points",
        "OS3 not counted (OS1 uses this amino-acid change)",
        "OS1_Moderate not counted (same change by c.1A>G: Likely oncogenic, VCV1; OS3 + OS2 for this change is "
        "stronger)",
        "OP4 not checked (absence from gnomAD not confirmed)",
        "OP4 not applied (absent from gnomAD, but few reads there)",
        "ambiguous (IntOGen: 12 mutated samples in 1 cohort with a role: too few to decide)",
        "OM4 not counted (OS1 uses this amino-acid change)",
        "OS1_Moderate not counted (same change by c.1A>G: Likely oncogenic, VCV1; OM3 for this change counts as much "
        "and is the standard's)",
        "OS2, OM4_Supporting not counted (OS1_Moderate uses this amino-acid change)",
        "OS1_Moderate not counted (same change by c.1A>G: Likely oncogenic, VCV1; OM4_Supporting + OP3 for this "
        "change is stronger)",
        "loss-of-function (IntOGen: 96% loss-of-function in 90 samples of cohorts with a role, 7 cohorts)",
        "ambiguous (IntOGen: 10 of 14 mutated samples in cohorts of ambiguous role)",
        "ambiguous (IntOGen: 38% loss-of-function in 13 samples, 5 cohorts; UniProt: tumour suppressor: they disagree)",
        "UniProt keywords unavailable: gene roles from IntOGen only",
        "Ensembl cannot be reached: only answers already in the cache are used",
        "PM2 met (absent from gnomAD, but few reads there: check)",
        "PM2 not assessed (no gnomAD data)",
        "NMD uncertain (the change spans an exon-intron boundary)",
        "PM2 not met (absent from gnomAD, but the panel requires 20x coverage)",
        "grpmax 2% is at or above its BA1 threshold (>= 0.4%): check its other BA1 conditions",
        "on ClinGen's BA1 exception list (Ghosh et al. 2018): BA1 does not apply, whatever the panel's threshold",
        "on the panel's BA1 exclusion list (well-known pathogenic variants): BA1 does not apply",
        "grpmax allele frequency 0.04% is above its BA1 threshold (>= 0.00486%), which it applies to the filtering "
        "allele frequency (not computed here): check the filtering allele frequency and its other BA1 conditions "
        "before BA1 (its other BA1 conditions: “x”)",
        "has no PM5 below Moderate: the tool lowered PM5 because the ClinVar reference is likely pathogenic or has one "
        "review star, and the panel judges the reference under its own specification",
        "the standard's rule for OS1 would give 6 points",
        "common: BA1 (grpmax 6.2%; allele numbers not checked in gnomAD)",
        "frequent: BS1 candidate (grpmax 7.04%); BA1 not met: the gene's expert panel sets BA1 above 10%",
        "OM4 not applied (this change is milder than the oncogenic ones at the residue)",
        "OM2 not applied (not a known cancer gene: UniProt, IntOGen)",
        "OM4_Supporting not used (p.Gly719Ser: Likely oncogenic; the hotspot criterion OM3 is stronger)",
        "not assessed (cancerhotspots.org or IntOGen unavailable)",
        "OS1 (same change by c.524_525delinsAT: Oncogenic, VCV000000001)",
        "OS2/SBS2 not checked (CIViC unavailable)",
        "OS2/SBS2 not applied (CIViC functional evidence conflicts)",
        "CIViC functional evidence rated below 3 stars not used",
        "not assessed automatically: functional studies beyond CIViC, OM1 functional domain, OP2, tumour type",
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
