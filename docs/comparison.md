# MissensExplorer compared with other tools

Two questions: which other software does something similar, and how good MissensExplorer's call is
next to the predictors those programs rely on. Checked in November 2026. The full tables are in
[benchmarks/results/tool_comparison.md](../benchmarks/results/tool_comparison.md), produced by
`python benchmarks/compare_tools.py`.

## Similar software

| Software | What it does | Missense call / PP3-BP4 | 3D structure | Access |
|---|---|---|---|---|
| [ProtVar](https://www.ebi.ac.uk/ProtVar) (EMBL-EBI) | Maps genomic or protein variants to UniProt and annotates them | Shows AlphaMissense, CADD, ESM1b, EVE, popEVE, conservation, FoldX, Missense3D; no combined call, no ACMG | Yes (PDB, AlphaFold, pockets, interfaces) | Free; REST API |
| [VarSome](https://varsome.com) | Variant search with automated ACMG classification | PP3/BP4 from MetaRNN | 3D viewer | Free tier 100 queries/day; API paid |
| [Franklin](https://franklin.genoox.com) (Genoox/QIAGEN) | ACMG interpretation platform | PP3/BP4 from REVEL | — | Free with registration; API in Premium |
| [GeneBe](https://genebe.net) | ACMG calculator (17 of 28 criteria) | PP3/BP4 from MetaRNN, else BayesDel | — | Free for research; API |
| [InterVar](https://wintervar.wglab.org) | Applies the 28 ACMG 2015 criteria, on ANNOVAR | PP3/BP4 from MetaSVM, GERP++, dbscSNV (not calibrated) | No | Source available, non-commercial |
| [OpenCRAVAT](https://opencravat.org) | Modular local annotator | Modules for AlphaMissense, REVEL, CADD, ESM1b, EVE, VARITY…; no combined call | Link-out | Open source (MIT) |
| [Ensembl VEP](https://www.ensembl.org/vep) | Consequence annotation with plugins | Shows SIFT, PolyPhen, AlphaMissense, CADD, REVEL, dbNSFP; no combined call | No | Open source; REST API |
| [Missense3D](https://missense3d.bc.ic.ac.uk) | Structural impact on PDB or AlphaFold models | Structural damaging/neutral flags | Yes | Free web; no API |
| [MutationTaster2021](https://www.genecascade.org/MutationTaster2021/) | Random-forest disease-causing prediction | Its own model; ClinVar/HGMD variants called by look-up | No | Free web |
| [Alamut Visual Plus](https://www.sophiagenetics.com/platform/alamut-visual-plus/) | Desktop interpretation software | ACMG points; REVEL, PolyPhen-2, splicing tools | — | Paid |

**Where MissensExplorer differs.**
- It combines AlphaMissense, REVEL, CADD and popEVE in one fitted model and shows each predictor's
  contribution to the call. None of the tools above has a combined model with that kind of breakdown.
- It shows the AlphaFold structure and the protein-level consequence, as ProtVar does.
- It gives PP3/BP4 from its own score, calibrated with ClinGen's method on ClinGen's data
  ([details](model.md#acmg-pp3bp4-from-the-model-score)), and keeps frequency as a separate, visible rule.
- It is open source, with a desktop app and a command line.
- It does **not** produce a full ACMG classification, as VarSome, Franklin, GeneBe or InterVar do.
  Those need evidence that a predictor cannot give: segregation, de novo status, functional studies.

**Why the platforms are not benchmarked directly.**
- VarSome and Franklin have no free bulk API, and their terms restrict automated use.
- Their missense evidence comes from one predictor each (MetaRNN, REVEL, MetaSVM, BayesDel), and
  those predictors are all in the benchmark below.

## Benchmark against other predictors

**Test set.** The held-out ClinVar benchmark: 1,972 expert-reviewed variants from 1,181 genes that
the MissensExplorer model never saw. 31 predictors were compared (33 scores, counting BayesDel with
and without AF and VARITY's leave-one-out version), their scores taken from dbNSFP through Ensembl.

**Measure.** AUC, the chance that a random pathogenic variant scores higher than a random benign
one. A difference is called real only when its 95% paired bootstrap interval excludes zero.

| AUC | All variants | Rare (< 0.1%) | New since 2023 |
|---|---|---|---|
| MetaRNN (ClinVar-trained, **uses frequency**) | **0.996** | **0.983** | **0.983** |
| ClinPred (ClinVar-trained, uses frequency) | 0.993 | 0.973 | 0.971 |
| BayesDel with AF (uses frequency) | 0.992 | 0.972 | 0.969 |
| MutationTaster2021 (ClinVar/HGMD look-up) | 0.992 | 0.974 | 0.972 |
| **MissensExplorer v2.5** | **0.985** | **0.974** | **0.971** |
| PHACTboost | 0.983 | 0.957 | 0.956 |
| REVEL (Franklin's PP3/BP4) | 0.982 | 0.972 | 0.966 |
| BayesDel without AF | 0.977 | 0.962 | 0.950 |
| gMVP | 0.972 | 0.966 | 0.956 |
| VARITY_R | 0.971 | 0.960 | 0.954 |
| AlphaMissense | 0.964 | 0.953 | 0.957 |
| MetaSVM (InterVar's PP3/BP4) | 0.954 | 0.916 | 0.907 |
| ESM1b | 0.935 | 0.920 | 0.909 |
| CADD | 0.923 | 0.885 | 0.890 |
| gnomAD frequency alone (reference) | 0.978 | 0.909 | 0.938 |

*New since 2023*: 454 variants whose pathogenic/benign label was not yet in the January 2023 ClinVar
release, so no tool trained before then can have seen it.

### What it shows

- **Best of the tools that leave out frequency.**
  - It significantly beats REVEL on all variants. On rare and new variants it is ahead by
    0.004–0.005, within noise.
  - It significantly beats AlphaMissense, BayesDel without AF, VARITY and PHACTboost on rare and new
    variants, and gMVP on new variants.
  - It also scores the 4% of variants that REVEL leaves out, through AlphaMissense.
- **MetaRNN ranks higher.**
  - MetaRNN is significantly ahead on all three sets (+0.009 to +0.012).
  - MetaRNN takes allele frequency as input. Adding frequency to the MissensExplorer score as a second
    feature gives the same AUC as MetaRNN: 0.996, 0.983 and 0.982. The tool offers that combination
    as an extra column for ranking variants
    ([score with frequency](model.md#score-with-frequency-for-ranking)). It is not used for the
    verdict or for ACMG evidence.
  - Whether its lead comes from frequency is not settled by these data. Among the 920 variants below
    0.01% or absent from gnomAD (865 pathogenic, 55 benign), MetaRNN's lead is not significant (0.961
    against 0.947; +0.013, 95% CI −0.008 to 0.042). Among the 412 present below 0.01% it is (0.966
    against 0.933; 95% CI of the difference 0.007 to 0.068), but frequency alone still ranks those
    variants with an AUC of 0.707, so the stratum does not remove frequency's information. MetaRNN was
    trained on ClinVar, which this benchmark samples; on laboratory data (below) it falls behind.
- **ClinPred, BayesDel with AF and MutationTaster2021** are ahead on all variants but level with
  MissensExplorer on rare and new ones. Their edge sits in the common benign variants, which
  frequency alone already separates (AUC 0.978).
- **Why MissensExplorer keeps frequency out of the score.** In ACMG/AMP classification, frequency has
  its own criteria (BA1, BS1, PM2). ClinGen therefore recommends computational evidence from "in
  silico tools that do not make direct use of allele frequency" ([Pejaver et al.
  2022](https://pmc.ncbi.nlm.nih.gov/articles/PMC9748256/)); otherwise the same frequency is counted
  twice. Frequency on ClinVar is also partly circular, since curators use it to call variants benign.
  The tool applies frequency as a separate rule; with it, the default verdict's AUC is 0.986.

### On laboratory data (independent of ClinVar)

Deep mutational scanning data from [ProteinGym](https://proteingym.org): 27,012 variants measured in 71
assays of 61 human proteins, each tool compared with MissensExplorer on the same variants
([full results](../benchmarks/results/proteingym_validation.md)). The order changes:

| Spearman with measured loss of function | | Difference from MissensExplorer (95% CI) |
|---|---|---|
| **MissensExplorer v2.5** | **0.462** | — |
| popEVE | 0.453 | −0.009 (−0.025 to +0.005) |
| AlphaMissense | 0.453 | −0.013 (−0.025 to −0.001) |
| MetaRNN | 0.431 | −0.034 (−0.052 to −0.018) |
| VARITY_R | 0.430 | −0.033 (−0.053 to −0.015) |
| REVEL | 0.419 | −0.043 (−0.053 to −0.032) |
| ClinPred | 0.402 | −0.059 (−0.073 to −0.044) |
| CADD | 0.337 | −0.125 (−0.152 to −0.100) |

- MissensExplorer is level with popEVE and slightly ahead of AlphaMissense, which use no clinical
  labels.
- MetaRNN and ClinPred, ahead on ClinVar, fall significantly behind it.

### Calls on rare variants

Each predictor used alone at its ClinGen-calibrated *supporting* thresholds, against the
MissensExplorer verdict:

| Rare variants | Pathogenic called damaging | Benign called damaging | Decisive calls | Decisive calls right |
|---|---|---|---|---|
| **MissensExplorer v2.5** | **91.9%** | 6.0% | **92.7%** | 96.9% |
| BayesDel without AF | 88.2% | 7.0% | 88.0% | 96.7% |
| REVEL | 86.0% | 6.0% | 85.3% | **97.6%** |
| VARITY_R | 77.7% | 5.0% | 77.4% | 96.7% |
| MutPred2 | 74.4% | 6.0% | 79.0% | 95.0% |
| ESM1b | 74.3% | 7.5% | 78.8% | 92.8% |
| AlphaMissense | 71.5% | **4.0%** | 74.5% | 95.3% |

- Against REVEL at the same false-alarm rate, MissensExplorer finds more pathogenic variants and
  gives a decisive call more often.
- REVEL calls fewer pathogenic variants benign (1.2% against 2.2%), so it is right a little more
  often when it decides.

## Caveats

- **All labels come from ClinVar.** Tools trained on ClinVar may have seen some of these variants;
  the *new since 2023* set limits that. HGMD and UniProt, used by others, overlap ClinVar's
  pathogenic variants.
- **The calibrations are older than some scores.** The ClinGen thresholds were set on older versions
  of some tools: CADD before v1.7, and ESM1b before dbNSFP updated its scores. Treat those rows as
  indicative.
- **dbNSFP carries the 2018 PrimateAI,** not PrimateAI-3D, so PrimateAI-3D is not tested here.
- **This benchmark was already used for MissensExplorer's design.** It showed the v2.2 problem, so it
  is not a fully fresh test. The ProteinGym check is: it was run once, after the design was fixed.
- **Small sets.** With 199 rare benign and 172 new benign variants, differences below about 0.01 AUC
  are within noise.
