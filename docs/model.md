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

Compared with 31 other predictors on the same variants, see [comparison.md](comparison.md).

The benchmark had already shown the v2.2 problem, so this is its second use for the model. The
benchmark took no part in choosing v2.3's design or cut-offs, which were set by the rule above on
the training genes alone.

## Independent validation on laboratory data (ProteinGym)

[Full results](../benchmarks/results/proteingym_validation.md), from `benchmarks/validate_proteingym.py`.

**Data.** Deep mutational scanning measures, in the laboratory, how each amino-acid substitution
changes a protein's function. The check used 13,408 single-nucleotide missense variants from 71
human assays on 61 proteins of [ProteinGym](https://proteingym.org). No ClinVar label is involved.

**Measure.** Spearman correlation with the measured loss of function, averaged per protein.

| | Spearman | vs MissensExplorer |
|---|---|---|
| **MissensExplorer v2.3** | **0.469** | — |
| popEVE | 0.460 | −0.009, not significant |
| AlphaMissense | 0.458 | −0.012, borderline |
| MetaRNN | 0.438 | −0.035, significant |
| REVEL | 0.425 | −0.044, significant |
| ClinPred | 0.405 | −0.065, significant |
| CADD | 0.344 | −0.126, significant |

- **The model generalizes beyond ClinVar.** It ranks first of the 34 scores tested.
- **It beats every input but one.** It does better than REVEL and CADD, and level with popEVE and
  AlphaMissense.
- **The ClinVar-trained meta-predictors drop back.** MetaRNN and ClinPred lead on ClinVar, where
  allele frequency helps them, but fall behind here.
- **Correlations are moderate for every tool,** because an assay measures one function in one
  experimental system, not disease.

## Somatic activating mutations: EGFR (outside the model's scope)

[Full results](../benchmarks/results/egfr_validation.md), from `benchmarks/validate_egfr.py`.

**Data.** Ikeuchi et al. (Oncogene 2022) tested every EGFR variant seen in at least two tumours in
COSMIC (most of them rare) for transformation of NIH-3T3 and Ba/F3 cells: 258 substitutions, plus
in-frame deletions and insertions. None is in the model's training set.

**Question.** The model separates disease-causing germline variants from harmless ones. Tumour
drivers such as EGFR's are *activating*, and a change that breaks the kinase is damaging without
activating it. How far does a damaging call carry over?

| Against | Variants | MissensExplorer AUC | Best other score |
|---|---|---|---|
| 3T3 transformation (score 3-4 vs 1) | 218 | 0.54 (0.43-0.64) | AlphaMissense 0.55 |
| OncoKB (oncogenic vs likely neutral) | 85 | 0.66 (0.49-0.81) | popEVE 0.68 |
| Ba/F3 transformation | 251 | 0.21 (0.15-0.28) | CADD 0.30 |

- **Activation is not what any of the scores measure.** Against 3T3 transformation every score is
  close to chance. The model calls 30 of 38 transforming variants damaging, and 127 of 187
  non-transforming ones too.
- **Loss of function is what they measure.** The 58 substitutions that failed in Ba/F3 also failed
  in 3T3. Most sit in the kinase domain, many on its catalytic motifs (glycine-rich loop, HRD loop,
  DFG). The model calls 54 of them damaging: these are broken kinases. Hence an AUC below 0.5.
- **The well-known uncommon lung-cancer mutations are flagged.** G719A/C/D/S, L861Q/R, S768I,
  E709A/K, L747S and T790M are called damaging, like L858R. L747P, a two-base change without
  AlphaMissense or REVEL scores, is damaging on a single predictor. Of the other variants the paper
  marks as uncommon, L833V, V769L/M and A864T are uncertain, and H870R is benign.
- **Use.** For a somatic variant, a damaging call means a change that matters to the protein, not an
  activating one. For oncogenicity, see the points below, with curated sources and functional data.

## Somatic oncogenicity points (ClinGen/CGC/VICC)

Optional (`--somatic`, or the *Somatic oncogenicity* switch). [Full results](../benchmarks/results/oncogenicity_validation.md),
from `benchmarks/validate_oncogenicity.py`; the rules are in `src/missensexplorer/oncogenicity.py`.

**What it is.** The ClinGen/CGC/VICC standard (Horak et al., Genet Med 2022) scores how likely a
somatic variant is to drive a tumour. Each criterion adds points: very strong 8, strong 4,
moderate 2, supporting 1, and the benign criteria subtract. The total gives the class:

| Points | 10 or more | 6 to 9 | 0 to 5 | -1 to -6 | -7 or fewer |
|---|---|---|---|---|---|
| Class | Oncogenic | Likely oncogenic | VUS | Likely benign | Benign |

These class limits are as restated by OncoVI (Carta et al., J Mol Diagn 2026). They match the two
worked examples in Horak et al.

**What is computed, from open data.**

| Criteria | Source |
|---|---|
| OS3 / OM3 / OP3: tumour hotspot, by the number of tumours at the residue and with this change | cancerhotspots.org |
| OS1 / OM4: the same change, or another missense change at the residue, classified oncogenic in ClinVar (one strength lower when only *likely* oncogenic) | ClinVar |
| OVS1: null variant in a loss-of-function gene; OM2: in-frame change in a cancer gene | IntOGen driver genes (CC0) |
| OP4: absent from gnomAD or below 1e-5; SBS1 / SBVS1: above 1% / 5% in a continental population | gnomAD |
| OS2 / SBS2: functional studies of this protein change (rated at least 3 of 5): evidence of an oncogenic effect, or of unaltered function | CIViC (CC0) |
| OP1 / SBP1: the model calls the missense variant damaging / benign | MissensExplorer |

For OS2, what counts as an oncogenic effect depends on the gene. Gain of function or neomorphic
function counts in any driver. Loss of function or a dominant-negative effect counts only in a
loss-of-function gene.

**Clinical evidence (no points).** CIViC's predictive, diagnostic and prognostic evidence for the
protein change goes in `clinical_evidence`. It gives the best level (A validated or in guidelines,
B clinical trials, C case reports, D preclinical, E inferential), the drug and the tumour type. This
is the other axis of somatic interpretation, what the variant means for treatment (AMP/ASCO/CAP),
so it never adds oncogenicity points. For example, L858R reaches level A (erlotinib in lung cancer)
and BRAF V600E level A (vemurafenib in melanoma), although CIViC has no functional study of either.

**What is not computed.**
- Functional studies outside CIViC, a critical functional domain (OM1) and OP2 need a curator. CIViC
  holds few functional studies: 341 in 63 genes, 196 of them rated 3 stars or more.
- The variant's own ClinVar oncogenicity classification never counts as evidence for itself; it is
  shown apart (`clinvar_oncogenicity`).
- OncoKB is not used: its terms do not allow use in clinical or commercial settings without a
  licence.

**Three choices.**
1. **OP1/SBP1 only outside activating genes.** In a gene IntOGen calls activating, the model's call
   does not count. The data:
   - In EGFR, an oncogene, a damaging call marks broken kinases as much as activated ones. With OP1
     the points ranked Ba/F3 transformation at 0.42, without it at 0.56.
   - In TP53, a tumour suppressor, the model is the strongest criterion: the points rank loss of
     function at 0.93 with it and 0.87 without it.
2. **OP4 below 1e-5.** The standard says "absent or extremely low" without a number.
3. **OS1 removes the hotspot criteria**, whose counts mix nucleotide changes. OM3 and OM4 are never
   combined; the stronger one is kept.

**Validation.** The points are ranked against laboratory assays: AUC, where ties count half.

| Gene (role) | Assay | Variants | Model alone | Points |
|---|---|---|---|---|
| EGFR (oncogene) | 3T3 transformation (3-4 vs 1), missense | 225 | 0.543 | **0.729** (0.620-0.823) |
| EGFR | Ba/F3 transformation, missense | 258 | 0.215 | **0.557** (0.507-0.608) |
| TP53 (tumour suppressor) | Kato yeast transactivation | 1,724 | 0.928 | **0.932** (0.917-0.946) |
| TP53 | Giacomelli loss of function | 2,158 | 0.841 | **0.842** (0.820-0.864) |

In this validation, CIViC items that cite the assays themselves, or yeast assays like Kato's, are left
out: 78 items. Otherwise the points would be tested against their own source. With the rest, OS2
applies to 3 EGFR and 10 TP53 variants.

- **Ranking:**
  - In the oncogene the points rank far better than the model, because the tumour hotspots carry
    the signal.
  - In the tumour suppressor they rank as well as the model, which is already good.
  - Without ClinVar the AUCs barely change: the points do not depend on ClinVar's classifications.
- **Calls:** they are cautious.
  - No TP53 variant functional in yeast is called oncogenic.
  - One EGFR variant that does not transform 3T3, K754E, comes out *Likely oncogenic*. CIViC holds
    evidence that it is oncogenic, and it does transform Ba/F3.
  - Few variants reach *Likely oncogenic* or *Oncogenic*:
    - TP53: 91 of 493 non-functional variants. G245S is *Oncogenic*, with 10 points from the
      hotspot, a CIViC functional study, the model and gnomAD.
    - EGFR: E746_A750del and T790M.
  - L858R, the TP53 hotspot R273H and most uncommon EGFR mutations stay *VUS*. They have no
    functional study in CIViC, and without one the automatic criteria rarely pass 5 points.
- **Use:** the points rank and document the evidence; they do not replace curation.
  - A *VUS* here means "not enough automatic evidence", not "probably harmless".
  - ClinVar's classification of the variant itself is shown next to them.
- **Not assessed:**
  - Insertions written as protein changes, such as most EGFR exon 20 insertions, cannot be traced
    to one DNA change; give the c. change.
  - Hotspots come from about 25,000 tumours (cancerhotspots.org), with less coverage of
    haematological, rare and paediatric cancers.

## Inputs that did not help

`benchmarks/feature_tests.py` adds one input at a time to the four predictors and measures the change
in AUC on rare training variants, with cross-validation split by gene.

**Rule.** Set before testing: an input is adopted only with a gain of at least 0.003, an interval
that excludes zero, and a coefficient of the expected sign.

**Result.** Nothing qualified:

| Added input | Gain in rare-variant AUC | Why it was rejected |
|---|---|---|
| Gene LoF constraint (gnomAD LOEUF) | +0.0030 | Wrong sign: it rewards LoF-tolerant genes |
| PrimateAI | +0.0024 | Below the rule's bar, and wrong sign next to the other predictors |
| Gene missense constraint (o/e, z) | +0.0006 to +0.0011 | Too small, and wrong sign |
| phyloP | +0.0006 | Too small, and wrong sign next to the other predictors |
| ESM1b, MutationAssessor | +0.0001 to +0.0002 | Too small |
| AlphaFold pLDDT, UniProt functional sites | +0.0003 | Too small (pLDDT available for half of the variants) |
| PROVEAN, SIFT4G | negative | — |

**Why gene constraint was rejected.**
- All three gene-constraint inputs took the wrong sign: they push towards *damaging* in genes that
  tolerate variation.
- In ClinVar, constrained genes carry many rare *benign* variants, found by clinical testing of
  dominant-disease genes. The model learns which genes get tested, not how damaging a substitution is.
- In the most constrained genes it would lower sensitivity on rare variants from 91.7% to 88.0%.
- Gene-level evidence also has its own ACMG criteria (PP2, BP1), so it would be counted twice.

**Earlier tests.** Non-linear models (splines, monotone boosting) and up-weighting rare benign variants
did not rank variants better either (see *Choice of cut-offs*). With the inputs available today, the
four-predictor logistic regression is at its ceiling.

## Limitations

- **ClinVar labels on both sides.** The model is trained and tested on ClinVar labels. REVEL was
  itself trained on disease-variant databases, and AlphaMissense's thresholds were calibrated on
  ClinVar. Absolute figures are therefore optimistic; the comparison between methods is the robust
  part. The check on laboratory data below does not depend on ClinVar.
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
- **ACMG PP3/BP4.** Computed from the model score, with thresholds calibrated by ClinGen's method;
  see the next section. `--acmg-source revel` (desktop app: *Options → PP3/BP4*) uses REVEL with the
  thresholds ClinGen published instead (AlphaMissense when REVEL is missing).

## ACMG PP3/BP4 from the model score

PP3 and BP4 are the ACMG/AMP criteria for computational evidence. ClinGen recommends taking them from
**one** score whose thresholds were calibrated with the method of
[Pejaver et al. (AJHG 2022)](https://pmc.ncbi.nlm.nih.gov/articles/PMC9748256/). That paper set the
thresholds of REVEL and other tools. MissensExplorer calibrates its own score with the same method,
on the same data, so its evidence has the same meaning as REVEL's.
`python benchmarks/calibrate_acmg.py` reproduces everything below; the full tables are in
[benchmarks/results/acmg_calibration.md](../benchmarks/results/acmg_calibration.md).

| Evidence (points) | Likelihood ratio needed | Model score |
|---|---|---|
| PP3 Strong (+4) | 33.5 | ≥ 0.985 |
| PP3 Moderate (+3) | 13.9 | ≥ 0.975 |
| PP3 Moderate (+2) | 5.8 | ≥ 0.941 |
| PP3 Supporting (+1) | 2.4 | ≥ 0.816 |
| BP4 Supporting (−1) | 1/2.4 | ≤ 0.271 |
| BP4 Moderate (−2) | 1/5.8 | ≤ 0.105 |
| BP4 Moderate (−3) | 1/13.9 | ≤ 0.017 |
| BP4 Strong (−4) | 1/33.5 | ≤ 0.005 |

The +3/−3 intervals count as Moderate under the 2015 guidelines, as for REVEL. The output reads, for
example, `PP3_Strong (model 0.994)` = +4 points.

**How the thresholds were set.**
- **Method.** A Python port of the authors' MATLAB code
  ([vpejaver/clingen-svi-comp_calibration](https://github.com/vpejaver/clingen-svi-comp_calibration)).
  - It assumes a prior probability of pathogenicity of 0.0441.
  - It estimates the posterior probability at each score in a window of at least 100 ClinVar
    variants and 3% of rare gnomAD variants.
  - The threshold of a level is where that posterior stops meeting the level's likelihood ratio.
    It is taken at the one-sided 95% bound over 10,000 bootstrap samples.
- **Check of the port.** On the paper's own data it reproduces the published thresholds of REVEL
  (0.644 / 0.774 / 0.879 / 0.932 for PP3; 0.290 / 0.181 / 0.052 / 0.016 for BP4), VEST4, MutPred2
  and PolyPhen-2. Every value is within 0.002, the size of bootstrap noise.
- **Data.**
  - The paper's ClinVar 2019 set: the data behind REVEL's ClinGen thresholds.
    - It has missense variants below 1% in gnomAD, with the training variants of REVEL and other
      tools removed.
    - 11,201 of its 11,834 variants are missense on the transcript MissensExplorer reports.
  - Each variant is scored by a model fitted without its gene, because the shipped model was trained
    on ClinVar variants of many of these genes. These scores agree with the shipped model's on genes
    neither model saw (correlation 0.9999).
  - The gnomAD part is 5,615 rare gnomAD v4 missense variants that are not in ClinVar.
- **Effect of the substitutions.** Done on REVEL's own scores:
  - **gnomAD set.** A set of 5,615 variants instead of 364,000 makes BP4 Moderate a little stricter
    (0.159 instead of 0.183) and leaves PP3 unchanged.
  - **Variant subset.** Keeping only the 11,201 variants that are missense on the current transcript
    moves BP4 Supporting from 0.290 to 0.342. The 633 variants left out include 344 pathogenic ones
    with middling REVEL scores, probably pathogenic for other reasons such as splicing.
    MissensExplorer gives PP3/BP4 only to variants that are missense on the reported transcript, so
    this is the right population. The tests below check the BP4 bands for exactly this reason.

**Tests on data not used for the calibration.** For each band, the likelihood ratio was compared
with the one the level needs:
- the paper's ClinVar 2020 set (8,723 variants), as in the paper;
- the rare variants of the held-out benchmark (1,163 variants, genes the shipped model never saw);
- among those, the ones labelled since 2023.

| ClinVar 2020 set | Pathogenic | Benign | Likelihood ratio (95% CI) | Needed |
|---|---|---|---|---|
| PP3 Strong | 1,233 | 50 | 55.6 (43.5–75.6) | 33.5 |
| PP3 +3 | 269 | 37 | 16.4 (12.0–23.8) | 13.9 |
| PP3 Moderate | 352 | 106 | 7.5 (6.1–9.4) | 5.8 |
| PP3 Supporting | 333 | 200 | 3.8 (3.2–4.5) | 2.4 |
| BP4 Supporting | 89 | 976 | 4.9 (4.0–6.1) | 2.4 |
| BP4 Moderate | 71 | 2,202 | 13.8 (11.1–17.9) | 5.8 |
| BP4 −3 | 10 | 982 | 43.5 (26.2–108.5) | 13.9 |
| BP4 Strong | 1 | 452 | 200 (62–∞) | 33.5 |

- **ClinVar 2020 set:** every band reaches the ratio it needs.
- **Held-out benchmark:** every band reaches it too.
- **Labelled since 2023:** one band falls short, PP3 +3 at 11.3 against 13.9. It rests on a single
  benign variant (31 pathogenic, 1 benign).

**Compared with REVEL at ClinGen's published thresholds** (the `revel` source):

| | ClinVar 2020: model | ClinVar 2020: REVEL | Held-out, rare: model | Held-out, rare: REVEL |
|---|---|---|---|---|
| Pathogenic with PP3 | **81.6%** | 73.2% | 87.8% | 88.5% |
| … PP3 Strong | **46.0%** | 28.2% | **46.5%** | 34.4% |
| Benign with BP4 | **76.3%** | 65.6% | **78.9%** | 71.4% |
| … BP4 Moderate or more | **60.2%** | 47.4% | **68.8%** | 51.3% |
| Pathogenic with BP4 (wrong way) | 6.4% | 6.7% | 2.2% | 1.5% |
| Benign with PP3 (wrong way) | 6.5% | 7.3% | 5.0% | 6.0% |

- The model gives stronger evidence more often, in both directions.
- Its wrong-way evidence is about as frequent as REVEL's:
  - ClinVar 2020 set: a little lower on both sides.
  - Held-out benchmark: 7 more pathogenic variants out of 964 get BP4 Supporting, and 2 fewer benign
    variants out of 199 get PP3.
- Among rare gnomAD variants, 20% get no computational evidence, against 28% with REVEL.

**Using it.**
- **Not published by ClinGen.** The calibration follows ClinGen's method on ClinGen's data, but it is
  this project's work, not reviewed by ClinGen. A laboratory that requires a ClinGen-published
  calibration should use `--acmg-source revel`.
- **One source only.** Use either the model or REVEL, chosen in advance, never both for the same
  variant. Pejaver et al. also recommend limiting PP3 and PM1 together to Strong (4 points).
- **Two different standards on one score.**
  - The verdict (*Likely damaging* from 0.76) is tuned to find pathogenic variants.
  - PP3 needs the likelihood ratio of the ACMG point system (from 0.816).
  - A variant between 0.76 and 0.816 is therefore *Likely damaging* without PP3. Likewise, BP4
    starts at 0.271 while *Likely benign* starts at 0.30.
- **Tied to one trained model.** The thresholds apply to model 2.3 as trained on 2026-11-23. After a
  retrain, the app falls back to REVEL until `benchmarks/calibrate_acmg.py` is run again, and a test
  fails as a reminder.

## PS1, PM5 and PM2

Three more ACMG/AMP criteria are computed automatically, next to PP3/BP4. They are listed with their
points in `acmg_codes`, and explained in `acmg_ps1_pm5` and `acmg_pm2`.

**PS1 and PM5: pathogenic variants at the same residue (ClinVar).**
- **Search.** For every missense SNV, the tool searches ClinVar (E-utilities, batched) for the
  pathogenic and likely pathogenic variants within two bases on each side, which covers the whole
  codon.
- **What counts as a reference variant.**
  - A germline classification of Pathogenic, Pathogenic/Likely pathogenic or Likely pathogenic, with
    review criteria (≥ 1 star, no conflicts).
  - The same MANE Select transcript, protein position and reference amino acid. ClinVar names its
    variants on that transcript.
  - Not the variant itself.
- **Strengths.**

  | Reference variant | Pathogenic | Only likely pathogenic |
  |---|---|---|
  | Same amino-acid change, other nucleotide change | PS1 (+4) | PS1_Moderate (+2) |
  | Another missense change at the residue | PM5 (+2) | PM5_Supporting (+1) |

  The lower strengths for likely pathogenic references follow several ClinGen expert panels (for
  example glaucoma and HNF4A). Only the stronger of PS1 and PM5 is given, never both.
- **Not given when SpliceAI ≥ 0.2.** The variant may then act through splicing rather than through
  the amino-acid change.
- **What you must still check.**
  - Whether the reference variant was classified independently of PS1/PM5.
  - Whether the reference variant itself acts through splicing (codons at exon edges).
  - Any gene-specific rule of a ClinGen expert panel; some, for example, give PM5_Strong when two
    pathogenic changes are known.
  - Codons split by an intron, or reference variants named on another transcript, are not found.

**PM2: absence from controls (gnomAD).**
- **Strength.** ClinGen lowered PM2 to Supporting in 2020 (SVI recommendation v1.0).
- **ClinGen's 2025 guidance on gnomAD v4.** Do not specify PM2 as strict absence; use a frequency
  threshold calculated for each disease, from its inheritance, prevalence and penetrance.
- **What the tool can decide.** It does not know the disease. But a variant that is truly absent
  meets every such threshold, so the tool gives **PM2_Supporting (+1) only for confirmed absence**:
  - gnomAD itself (v4 on GRCh38, v2.1 on GRCh37) does not have the variant;
  - and at least half of gnomAD's exome samples have ≥ 20 reads at the position.
- **Other cases.** A poorly covered position gives *not applicable*. A variant present in gnomAD
  gives *depends on the disease*: compare its frequency with the disease's maximum credible allele
  frequency.
- **Why the tool asks gnomAD directly.** Ensembl attaches gnomAD frequencies only to variants in its
  variation database, so many rare gnomAD variants, singletons above all, came back with no
  frequency.
  - Before this check, `population_frequency` called them *not found in gnomAD*.
  - The tool now asks gnomAD about every variant without a frequency, and fills in its real one.
  - In a random sample of 300 rare gnomAD v4 missense variants, 46% were in that situation: exactly
    the ones without a dbSNP identifier.
- **Combining with PP3/BP4.** Pejaver et al. (2022) allow PM2 and BS1 to be combined with PP3/BP4
  without limits, because the PP3/BP4 score does not use allele frequency.

**How informative they are** ([results](../benchmarks/results/acmg_criteria.md);
`python benchmarks/acmg_criteria.py`).
- **Data.** The held-out benchmark: 1,972 variants from genes the model never saw.
- **Time split.** The reference variants come from the January 2023 ClinVar release, and the test
  variants were classified after it.

| Criterion | Pathogenic with it | Benign with it | Likelihood ratio | Needed |
|---|---|---|---|---|
| Any PS1/PM5, current ClinVar | 50.9% | 0.3% | 169 | ≥ 2.4 |
| Any PS1/PM5, time split | 26.6% | 1.2% | 23 | ≥ 2.4 |
| PM5_Supporting, time split | 11.3% | 1.2% | 9.8 (95% CI from 3.5) | 2.4 |
| PM2_Supporting, rare variants | 46.7% | 1.0% | 46 (95% CI from 17) | 2.4 |

- Every criterion clears, by far, the likelihood ratio its strength needs.
- The ratios are optimistic. ClinVar submitters use these same criteria, and frequency, to classify
  variants, so the labels are not independent of the evidence. The time split limits this for
  PS1/PM5 but does not remove it.
- The strengths are therefore kept at the published ones, not raised.

## Score with frequency (for ranking)

`score_with_frequency` combines the model score with the gnomAD grpmax frequency into one number for
ranking the variants of a sample. Desktop app: column *With frequency*.

```
z = −6.937 + 0.855·logit(model score) − 1.822·log10(gnomAD grpmax + 10⁻⁶)
score with frequency = 1 / (1 + e^(−z))
```

**How it was made.**
- **Specification.** Fixed earlier, as the exploratory *frequency as a feature* row of the tool
  comparison. It is a logistic regression on the model's training set (9,218 variants).
- **Not retuned.** So the held-out benchmark is not used a second time to choose it.
- **Absent variants.** A variant absent from gnomAD's main ancestry groups counts as frequency 0,
  that is about one allele in a million.
- **Reproduce.** `python benchmarks/train_frequency_score.py`; the
  [results](../benchmarks/results/frequency_score.md) are in the repository.

| AUC, held-out benchmark | Model alone | Score with frequency | Frequency alone | MetaRNN |
|---|---|---|---|---|
| All variants (1,969) | 0.985 | **0.996** (+0.010; 95% CI +0.007 to +0.014) | 0.978 | 0.996 |
| Rare, < 0.1% (1,163) | 0.974 | **0.983** (+0.009; +0.003 to +0.014) | 0.909 | 0.983 |
| Labelled since 2023 (454) | 0.971 | **0.982** (+0.011; +0.003 to +0.022) | 0.938 | 0.983 |

With frequency, the score separates pathogenic from benign variants as well as MetaRNN, the best
tool in the comparison.

**Read it as a ranking, not as a probability or a call.**
- **Frequency weighs heavily, the way ClinVar's variants do.** A variant absent from gnomAD gets a
  high value even with a low model score: 0.98 at model score 0.5, 0.89 at 0.1, and about 0.5 at
  0.01. In ClinVar most absent variants are pathogenic. In a patient's exome, most absent variants
  are private and benign.
- **Use it to sort a sample's variants.** It brings forward the rare and absent variants with high
  model scores. Do not compare it with the model's cut-offs (0.76 / 0.30), and do not cut it at a
  threshold.
- **Not ACMG evidence.** Frequency has its own criteria (BA1, BS1, PM2), and ClinGen asks for PP3/BP4
  from a score that leaves frequency out, which the model score does. Using this score for PP3/BP4
  would count frequency twice.
- **Common pathogenic variants rank low.** Variants that are frequent but pathogenic (recessive or
  founder conditions, for example HFE C282Y) rank low here. The main prediction keeps them, through
  ClinGen's BA1 exception list.
- **Tied to one trained model.** The coefficients belong to model 2.3 as trained on 2026-11-23. After
  a retrain the column stays empty until `train_frequency_score.py` is run again, and a test fails
  as a reminder.
