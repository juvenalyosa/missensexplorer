# Missense model (version 2.3)

MissensExplorer calls missense variants with a logistic regression over four predictor scores. This
page describes what the model is, how it was trained and validated, and where it should not be
trusted. Retrain it with `python benchmarks/train_model.py` (needs `pip install -e ".[benchmark]"`).
The coefficients live in `src/missensexplorer/model_params.py`.

## What it computes

```
z = 0.0255 + 1.161·AlphaMissense′ + 1.661·REVEL′ + 0.897·CADD′ − 0.576·popEVE′
model score = 1 / (1 + e^(−z))
```

Each `′` is the score standardised with the training mean and standard deviation. A missing score
is replaced by the training median, so it does not push the call either way. The model is used only
when AlphaMissense or REVEL has a score. Otherwise the v2.1 weighted vote decides.

| Model score | Label |
|---|---|
| ≥ 0.76 | *Likely damaging* |
| ≤ 0.30 | *Likely benign* |
| in between | *Uncertain (conflicting predictors)* when the predictors disagree, otherwise *Uncertain (intermediate scores)* |

**The score is a ranking, not a probability.** It was learned with the two classes weighted
equally, as REVEL was. It is *not* the chance that a variant is pathogenic in a given patient, which
also depends on the gene, the disease and how the variant was found.

**Every call can be explained.** Each predictor's contribution is its term in `z`, in log-odds;
positive pushes towards damaging. The contributions appear in `prediction_basis` (for example
`model score 0.99 (REVEL +2.3, AlphaMissense +1.5, CADD +0.7)`) and as a chart in the desktop app.

## Training data

- **Source:** ClinVar `variant_summary.txt.gz`, as in the benchmark.
- **Selection:** germline missense SNVs on GRCh38 with ≥ 2 review stars. Pathogenic/Likely pathogenic
  versus Benign/Likely benign. At most 8 variants per gene and class, so a few genes cannot dominate.
- **Held out:** every gene of the benchmark sample (1,198 genes) was excluded, so the benchmark below
  measures genes that the model never saw.
- **Enriched in rare benign variants:** 9,218 variants from 3,263 genes:
  - every pathogenic variant (4,609);
  - every *rare* benign variant, with gnomAD grpmax < 0.1% or absent (1,826);
  - 2,783 common benign variants, picked at random to balance the classes.
- **Annotation:** Ensembl VEP REST, release 116, with the same pipeline that the tool runs.

**Why enriched.** Rare benign variants are the hard cases, and the ones the model is used for:
common variants are handled by the frequency rule, and most of the rare missense variants in a
patient are benign. In a plain random ClinVar sample only 12% of the benign variants are rare,
against 98% of the pathogenic ones, so a model trained on it learns "benign" mostly from common
variants. Version 2.2 was trained that way, and on rare variants it called more benign variants
damaging than the v2.1 vote did (benchmark: 10.6% against 9.0%).

The training script refuses to fit when any annotation source failed, or when fewer than 90% of
the variants have an AlphaMissense or REVEL score. A first v2.2 run that lost the connection halfway
gave a useless model; this check exists because of it.

## Choice of cut-offs

The cut-offs are chosen by the training script with 5-fold cross-validation split by gene, on the
6,355 rare variants of the training set (the benchmark plays no part):

- **Damaging: 0.76.** The highest cut-off (at least 0.70) at which the model still calls at least
  as many rare pathogenic variants damaging as the v2.1 vote did. Within that limit it calls as few
  rare benign variants damaging as possible, because false alarms are the costlier error among a
  patient's rare variants. Variants just below the cut-off become *uncertain*, never *benign*.
- **Benign: 0.30.** The highest cut-off (at most 0.30) at which no more rare pathogenic variants
  are called benign than with the v2.1 vote.

Cross-validated on the rare training variants:

| | v2.1 vote | model at 0.70/0.30 | **model at 0.76/0.30** |
|---|---|---|---|
| Pathogenic called *damaging* | 84.9% | 87.1% | **84.9%** |
| Benign called *damaging* | 7.1% | 6.9% | **5.8%** |
| Pathogenic called *benign* | 5.5% | 4.9% | **4.9%** |

The same cross-validation compared other designs; none ranked rare variants better, so the
simplest was kept:

| Training design (rare variants, CV by gene) | AUC |
|---|---|
| v2.2: random sample, 360 rare benign variants | 0.964 |
| **v2.3: enriched, 1,826 rare benign variants** | **0.965** |
| Rare variants only | 0.964 |
| Additive splines (non-linear per predictor) | 0.964 |
| Gradient boosting, monotone | 0.965 |
| v2.1 vote fraction | 0.933 |

Up-weighting the rare benign variants further only moves the operating point (the same as a
higher cut-off); it does not rank them better.

## Choice of predictors

All seven available scores were tried on the v2.2 training set, using 5-fold cross-validation split
by gene:

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
genes, none used for training.

| All variants | v2.1 vote | v2.2 model | **v2.3 model** |
|---|---|---|---|
| AUC | 0.959 | 0.985 | **0.985** |
| Pathogenic called *damaging* | 89.1% | 93.1% | **91.4%** |
| Benign called *benign* | 87.3% | 88.4% | **88.4%** |
| Pathogenic called *benign* | 2.7% | 2.4% | **2.4%** |
| Benign called *damaging* | 5.2% | 5.1% | **3.2%** |
| Decisive calls | 92.1% | 94.5% | **92.7%** |
| Decisive calls that are right | 95.7% | 96.0% | **97.0%** |

The v2.2 and v2.3 columns include the frequency rule (BA1); the v2.1 column is the vote alone.

| Rare variants (gnomAD grpmax < 0.1%: 964 pathogenic, 199 benign) | v2.1 vote | v2.2 model | **v2.3 model** |
|---|---|---|---|
| AUC | 0.937 | 0.973 | **0.974** |
| Pathogenic called *damaging* | 89.2% | 93.4% | **91.6%** |
| Benign called *damaging* | 9.0% (18) | 10.6% (21) | **6.0% (12)** |
| Pathogenic called *benign* | 2.6% | 2.3% | **2.3%** |
| Benign called *benign* | 81.9% | 80.4% | **79.4%** |
| Decisive calls that are right | 96.0% | 96.1% | **96.8%** |

- On rare variants, v2.3 makes fewer wrong calls than v2.1 in both directions and finds more of the
  pathogenic variants.
- The price is that more rare benign variants stay *uncertain* (14.6% against 9.1% with v2.1), so
  slightly fewer are called *benign*.
- Compared with v2.2, the borderline variants (scores 0.70–0.76) are now *uncertain* instead of
  *damaging*. Sensitivity drops by about two points, and about 40% fewer benign variants are
  called damaging.
- **Variants evaluated since 2024:** AUC 0.984 (v2.1: 0.960). 3.4% of the benign variants are
  called damaging (v2.1: 5.1%).

The benchmark had already shown the v2.2 problem, so this is its second use for the model. The
benchmark took no part in choosing v2.3's design or cut-offs, which were set by the rule above on
the training genes alone.

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
  - With a 1% cut-off, specificity in the benchmark rises to 91.4%. Part of that gain is circular,
    because curators use the same frequencies to call variants benign. In the 5,874-variant run,
    13 pathogenic variants were ≥ 1% (CYP21A2, HBB E7K, TTR V142I, BBS1 M390R …), mostly in
    recessive conditions. That is why the default is 5%.
- **ACMG PP3/BP4.**
  - Uses the ClinGen-calibrated intervals of one predictor: REVEL (Pejaver et al., AJHG 2022), or
    AlphaMissense when REVEL is missing (Bergquist et al., Genet Med 2025).
  - Points: Supporting ±1, Moderate ±2 (±3 for the intermediate interval), Strong ±4.
  - Used alone, it is right more often when it decides (97.2%) but decides less often (84.8%) than
    the model.
