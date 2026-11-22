# Missense model (version 2.2)

MissensExplorer calls missense variants with a logistic regression over four predictor scores. This
page describes what the model is, how it was trained and validated, and where it should not be
trusted. Retrain it with `python benchmarks/train_model.py` (needs `pip install -e ".[benchmark]"`).
The coefficients live in `src/missensexplorer/model_params.py`.

## What it computes

```
z = 0.0107 + 1.150·AlphaMissense′ + 1.752·REVEL′ + 0.838·CADD′ − 0.828·popEVE′
model score = 1 / (1 + e^(−z))
```

Each `′` is the score standardised with the training mean and standard deviation. A missing score
is replaced by the training median, so it does not push the call either way. The model is used only
when AlphaMissense or REVEL has a score. Otherwise the v2.1 weighted vote decides.

| Model score | Label |
|---|---|
| ≥ 0.70 | *Likely damaging* |
| ≤ 0.30 | *Likely benign* |
| in between | *Uncertain (conflicting predictors)* when the predictors disagree, otherwise *Uncertain (intermediate scores)* |

**The score is a ranking, not a probability.** It was learned from equal numbers of pathogenic and
benign variants, as REVEL was. It is *not* the chance that a variant is pathogenic in a given
patient, which also depends on the gene, the disease and how the variant was found.

**Every call can be explained.** Each predictor's contribution is its term in `z`, in log-odds;
positive pushes towards damaging. The contributions appear in `prediction_basis` (for example
`model score 0.97 (REVEL +2.1, AlphaMissense +1.5, CADD +0.8)`) and as a chart in the desktop app.

## Training data

- **Source:** ClinVar `variant_summary.txt.gz`, as in the benchmark.
- **Selection:** germline missense SNVs on GRCh38 with ≥ 2 review stars. Pathogenic/Likely pathogenic
  versus Benign/Likely benign.
- **Size:** 5,837 variants from 2,682 genes, balanced: 2,851 pathogenic and 2,986 benign. At most
  8 variants per gene and class, so a few genes cannot dominate.
- **Held out:** every gene of the benchmark sample (1,198 genes) was excluded, so the benchmark below
  measures genes that the model never saw.
- **Annotation:** Ensembl VEP REST, release 116, with the same pipeline that the tool runs.

The training script refuses to fit when any annotation source failed, or when fewer than 90% of
the variants have an AlphaMissense or REVEL score. A first run that lost the connection halfway
gave a useless model; this check exists because of it.

## Choice of predictors

All seven available scores were tried, using 5-fold cross-validation on the training set split by
gene:

| Predictors | AUC |
|---|---|
| AlphaMissense, REVEL, CADD, popEVE, EVE, SIFT, PolyPhen-2 | 0.9766 |
| **AlphaMissense, REVEL, CADD, popEVE** | **0.9762** |
| AlphaMissense, REVEL, CADD | 0.9750 |

EVE, SIFT and PolyPhen-2 add nothing on top of the other four. Because they overlap with them,
they also get coefficients of the "wrong" sign: a more tolerated SIFT score pushed towards
damaging. That would make the explanation misleading, so they were left out. They still appear in
the table, and SIFT/PolyPhen-2 still decide when nothing else has a score (labelled *limited
evidence*).

Other options did not help in an exploratory run on 5,874 variants:
- Gradient boosting had the same AUC as the logistic regression.
- Grantham, BLOSUM62 and SpliceAI added nothing.
- Missing-value indicators added nothing.

## Validation (held-out genes)

[Full benchmark](../benchmarks/results/clinvar_benchmark.md): 1,972 ClinVar variants from 1,181
genes, none used for training. The cut-offs were fixed at 0.7/0.3 before the benchmark was run.

| | v2.1 vote | v2.2 model | v2.2 default (model + BA1) |
|---|---|---|---|
| AUC | 0.959 | 0.985 | — |
| Pathogenic called *damaging* | 89.1% | 93.1% | 93.1% |
| Benign called *benign* | 87.3% | 88.0% | 88.4% |
| Pathogenic called *benign* | 2.7% | 2.4% | 2.4% |
| Benign called *damaging* | 5.2% | 5.1% | 5.1% |
| Decisive calls | 92.1% | 94.3% | 94.5% |
| Decisive calls that are right | 95.7% | 96.0% | 96.0% |

- **Rare variants** (gnomAD grpmax < 0.1%: 964 pathogenic, 199 benign):
  - the AUC rises from 0.937 to 0.973 and sensitivity from 89.2% to 93.4%;
  - 10.6% of the rare benign variants are called damaging, against 9.0% with v2.1. That is 21
    against 18 variants, a gap within sampling noise but not in the model's favour.
- **Variants evaluated since 2024:** AUC 0.983 (v2.1: 0.960).
- **Cross-validated on the training set:** AUC 0.976.

## Limitations

- **ClinVar labels on both sides.** The model is trained and tested on ClinVar labels. REVEL was
  itself trained on disease-variant databases, and AlphaMissense's thresholds were calibrated on
  ClinVar. Absolute figures are therefore optimistic; the comparison between methods is the robust
  part. An independent check against deep mutational scanning data (for example MaveDB) has not
  been done yet.
- **Well-studied genes are over-represented.** ClinVar variants with ≥ 2 stars cluster in such genes.
- **Germline missense only.** The model is not trained for somatic or cancer variants, or for other
  consequences; those keep their rules (loss of function, SpliceAI, CADD for in-frame changes).
- **A predictor update needs a retrain.** If Ensembl changes a predictor's version, retrain and
  rerun the benchmark. The run record of each analysis stores the model version and training date.

## Separate rules shown next to the model

- **Population frequency (BA1).**
  - A damaging or uncertain call becomes *Likely benign (common in population)* at gnomAD
    grpmax ≥ 5%. grpmax is the highest frequency among the African/African American, Admixed
    American, East Asian, non-Finnish European and South Asian groups.
  - The protein-level call stays in `prediction_basis`.
  - The nine variants on ClinGen's BA1 exception list (Ghosh et al., Hum Mutat 2018) are never
    called benign for being common. They are matched by dbSNP id and allele, so the transcript
    does not matter; BTD p.Asp444His, for example, is p.Asp424His on the MANE transcript. Without
    the list, HFE C282Y (grpmax 7.1%) would be called benign.
  - With a 1% cut-off, specificity in the benchmark rises to 91.5%. Part of that gain is circular,
    because curators use the same frequencies to call variants benign. In the 5,874-variant run,
    13 pathogenic variants were ≥ 1% (CYP21A2, HBB E7K, TTR V142I, BBS1 M390R …), mostly in
    recessive conditions. That is why the default is 5%.
- **ACMG PP3/BP4.**
  - Uses the ClinGen-calibrated intervals of one predictor: REVEL (Pejaver et al., AJHG 2022), or
    AlphaMissense when REVEL is missing (Bergquist et al., Genet Med 2025).
  - Points: Supporting ±1, Moderate ±2 (±3 for the intermediate interval), Strong ±4.
  - Used alone, it is right more often when it decides (97.2%) but decides less often (84.8%) than
    the model.
