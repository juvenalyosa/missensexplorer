# EGFR: predictions against a transformation assay

Data: Ikeuchi et al., Oncogene 2022, doi:10.1038/s41388-022-02263-4 (Supplementary Table 4), CC BY 4.0. Every EGFR variant seen in at least two tumours in COSMIC v81, tested for
transformation of NIH-3T3 cells (focus formation: 1 = none, 2 = mild, 3 = moderate, 4 = severe) and of
Ba/F3 cells (growth without IL-3). The table lists 356 variants and 39 with 3T3 score 2, where the paper's
text gives 357 and 29; the table is used as it is. Produced by `benchmarks/validate_egfr.py` on 2026-11-26 (Ensembl 116).

The model scores how likely a missense change is to cause *germline disease*; the EGFR variants that drive
tumours are *activating*. A variant that breaks the kinase is damaging but does not transform, so this
tests how far a damaging call carries over to an activating one. None of these variants is in the model's
training set (overlap: none).

## Missense variants

258 substitutions; 251 received a model score. Transforming = 3T3 score 3-4; not transforming = score 1; score 2 (mild) is left out of the AUCs.

| Score | Variants | Transforming | AUC (95% CI) | Difference from the model |
|---|---:|---:|---|---|
| MissensExplorer v2.5 | 218 | 34 | 0.542 (0.433-0.640) |  |
| AlphaMissense | 218 | 34 | 0.545 (0.439-0.645) | -0.040 to +0.050 |
| REVEL | 218 | 34 | 0.528 (0.422-0.639) | -0.060 to +0.028 |
| CADD | 218 | 34 | 0.500 (0.389-0.606) | -0.133 to +0.049 |
| popEVE | 218 | 34 | 0.537 (0.439-0.634) | -0.047 to +0.041 |

Inside the kinase domain (residues 712-979) only:

| Score | Variants | Positives | AUC (95% CI) |
|---|---:|---:|---|
| MissensExplorer v2.5 | 151 | 20 | 0.569 (0.433-0.694) |
| AlphaMissense | 151 | 20 | 0.612 (0.485-0.730) |
| REVEL | 151 | 20 | 0.546 (0.401-0.696) |
| CADD | 151 | 20 | 0.488 (0.366-0.610) |
| popEVE | 151 | 20 | 0.550 (0.432-0.674) |

Outside the kinase domain only:

| Score | Variants | Positives | AUC (95% CI) |
|---|---:|---:|---|
| MissensExplorer v2.5 | 67 | 14 | 0.567 (0.374-0.757) |
| AlphaMissense | 67 | 14 | 0.505 (0.333-0.681) |
| REVEL | 67 | 14 | 0.596 (0.424-0.758) |
| CADD | 67 | 14 | 0.564 (0.377-0.751) |
| popEVE | 67 | 14 | 0.599 (0.410-0.779) |

Against the Ba/F3 assay (transformed or not), all scores 1-4:

| Score | Variants | Positives | AUC (95% CI) |
|---|---:|---:|---|
| MissensExplorer v2.5 | 251 | 194 | 0.214 (0.154-0.280) |
| AlphaMissense | 251 | 194 | 0.217 (0.156-0.288) |
| REVEL | 251 | 194 | 0.248 (0.179-0.317) |
| CADD | 251 | 194 | 0.300 (0.230-0.374) |
| popEVE | 251 | 194 | 0.205 (0.143-0.275) |

Against OncoKB (Oncogenic or Likely Oncogenic vs Likely Neutral):

| Score | Variants | Positives | AUC (95% CI) |
|---|---:|---:|---|
| MissensExplorer v2.5 | 85 | 72 | 0.656 (0.490-0.811) |
| AlphaMissense | 85 | 72 | 0.641 (0.487-0.785) |
| REVEL | 85 | 72 | 0.665 (0.509-0.806) |
| CADD | 85 | 72 | 0.640 (0.457-0.818) |
| popEVE | 85 | 72 | 0.676 (0.528-0.815) |

Calls of the pipeline per 3T3 score (one row per protein change; *mixed* = its nucleotide changes got
different calls):

| 3T3 score | Damaging | Uncertain | Benign / silent | Not assessed | mixed | error |
|---|---:|---:|---:|---:|---:|---:|
| 4 | 12 | 0 | 3 | 0 | 0 | 0 |
| 3 | 18 | 3 | 2 | 0 | 0 | 0 |
| 2 | 18 | 6 | 9 | 0 | 0 | 0 |
| 1 | 127 | 28 | 30 | 0 | 2 | 0 |

- Transforming (3-4) called damaging: 30 of 38 (79%); called benign: 5 of 38 (13%).
- Not transforming (1) called damaging: 127 of 187 (68%); called benign: 30 of 187 (16%).

## Reading

- **Activation is not what the scores measure.** Against 3T3 transformation every score is close to chance
  (AUC about 0.5 to 0.55). A damaging call covers most transforming variants, but also most of the
  variants that do not transform.
- **Loss of function is.** The 58 substitutions that did not transform Ba/F3 cells all also failed in 3T3 (58 of 58 with score 1), 51 lie in the kinase domain, and 54 are called damaging (median model score 0.98, against 0.86 for the variants that transformed Ba/F3). Several sit on the kinase's catalytic motifs (glycine-rich loop GSGAFG (719-724): G721A, G721S; catalytic loop HRDLAARN (835-842): A839T, D837N, L838P, N842D, R836H, R836S, R841K; DFG motif (855-857): D855N, F856L, G857E): changes that break the kinase. Hence the Ba/F3 AUCs below 0.5: there, "does not transform" mostly means "broken", which the predictors detect.
- So, for a somatic EGFR variant, a damaging call says the protein is changed where it matters, not that it
  is activating; whether it drives a tumour needs curated or functional data (see also `--somatic`).
  A benign call does not rule activation out either: 5 of 38 (13%) of the transforming variants get one, against 30 of 187 (16%) of the others.

## The uncommon mutations of lung cancer

The variants the paper's table marks as *uncommon EGFR mutations* (the column also marks L858R and T790M).
3T3: focus formation score; Ba/F3: transformed (yes/no).

| Variant | COSMIC | 3T3 | Ba/F3 | DNA change | Call | Model | AlphaMissense | REVEL | CADD | popEVE |
|---|---:|---:|---|---|---|---:|---:|---:|---:|---:|
| L62R | 24 | 2 | yes | c.185T>G | Likely damaging | 0.93 | 0.74 | 0.80 | 26.2 | -3.46 |
| A702S | 2 | 1 | yes | c.2104G>T | Uncertain (conflicting predictors) | 0.51 | 0.37 | 0.49 | 25.4 | -3.60 |
| I706T | 4 | 3 | yes | c.2117T>C | Likely damaging | 0.98 | 0.96 | 0.85 | 26.5 | -3.84 |
| E709_T710-D | 33 | 3 | yes |  | not annotated (insertion in protein notation) |  |  |  |  |  |
| E709A | 37 | 3 | yes | c.2126A>C | Likely damaging | 0.93 | 0.90 | 0.59 | 27.2 | -3.91 |
| E709G | 18 | 2 | yes | c.2126A>G | Likely damaging | 0.96 | 0.87 | 0.67 | 29.0 | -4.00 |
| E709H | 2 | 3 | yes | c.2125_2127delinsCAT, c.2125_2127delinsCAC | Likely damaging (single predictor) |  |  |  |  |  |
| E709K | 63 | 4 | yes | c.2125G>A | Likely damaging | 0.92 | 0.96 | 0.47 | 29.1 | -4.04 |
| L718V | 2 | 1 | yes | c.2152C>G | Likely damaging | 0.80 | 0.79 | 0.54 | 21.9 | -3.88 |
| G719A | 206 | 3 | yes | c.2156G>C | Likely damaging | 0.98 | 0.99 | 0.82 | 25.9 | -4.04 |
| G719C | 93 | 4 | yes | c.2155G>T | Likely damaging | 0.99 | 1.00 | 0.88 | 28.9 | -4.30 |
| G719D | 16 | 3 | yes | c.2156G>A | Likely damaging | 0.99 | 1.00 | 0.86 | 26.9 | -4.24 |
| G719S | 159 | 3 | yes | c.2155G>A | Likely damaging | 0.99 | 1.00 | 0.85 | 27.6 | -4.19 |
| S720F | 16 | 2 | yes | c.2159C>T | Likely damaging | 0.85 | 0.66 | 0.60 | 25.5 | -3.97 |
| G721S | 3 | 1 | no | c.2161G>A | Likely damaging | 0.99 | 1.00 | 0.86 | 27.6 | -4.23 |
| G724S | 15 | 1 | yes | c.2170G>A | Likely damaging | 0.97 | 1.00 | 0.74 | 27.6 | -3.98 |
| T725M | 10 | 1 | yes | c.2174C>T | Likely damaging | 0.80 | 0.40 | 0.68 | 26.5 | -3.83 |
| L747_P753-Q | 0 | 4 | yes |  | not annotated (insertion in protein notation) |  |  |  |  |  |
| L747P | 27 | 4 | yes | c.2239_2240delinsCC | Likely damaging (single predictor) |  |  |  |  |  |
| L747S | 16 | 4 | yes | c.2240T>C | Likely damaging | 0.99 | 1.00 | 0.84 | 27.3 | -4.20 |
| P753S | 18 | 1 | yes | c.2257C>T | Uncertain (intermediate scores) | 0.49 | 0.39 | 0.51 | 24.9 | -3.35 |
| K754E | 6 | 1 | yes | c.2260A>G | Likely damaging | 0.89 | 0.74 | 0.64 | 26.6 | -3.81 |
| V765M | 11 | 1 | yes | c.2293G>A | Uncertain (conflicting predictors) | 0.65 | 0.68 | 0.29 | 32.0 | -3.50 |
| S768I | 273 | 4 | yes | c.2303G>T | Likely damaging | 0.98 | 0.90 | 0.77 | 31.0 | -3.96 |
| V769_D770insASV | 58 | 4 | yes |  | not annotated (insertion in protein notation) |  |  |  |  |  |
| V769_D770insGVV | 3 | 3 | yes |  | not annotated (insertion in protein notation) |  |  |  |  |  |
| V769L | 17 | 2 | yes | c.2305G>C, c.2305G>T | Uncertain (conflicting predictors) | 0.59 | 0.86 | 0.29 | 25.6 | -3.44 |
| V769M | 12 | 2 | yes | c.2305G>A | Uncertain (conflicting predictors) | 0.61 | 0.57 | 0.32 | 33.0 | -3.40 |
| N771-GY | 0 | 4 | yes |  | not annotated (insertion in protein notation) |  |  |  |  |  |
| H773_V774insH | 8 | 3 | yes |  | not annotated (insertion in protein notation) |  |  |  |  |  |
| H773L | 11 | 2 | yes | c.2318A>T | Likely damaging | 0.88 | 0.53 | 0.66 | 28.4 | -4.08 |
| V774M | 23 | 3 | yes | c.2320G>A | Likely damaging | 0.97 | 0.93 | 0.62 | 33.0 | -3.96 |
| R776C | 9 | 1 | yes | c.2326C>T | Likely damaging | 0.96 | 0.84 | 0.59 | 34.0 | -3.92 |
| R776H | 21 | 1 | yes | c.2327G>A | Likely damaging | 0.95 | 0.65 | 0.68 | 34.0 | -3.81 |
| V786M | 12 | 1 | yes | c.2356G>A | Uncertain (intermediate scores) | 0.70 | 0.79 | 0.37 | 27.3 | -3.61 |
| T790M | 1398 | 3 | yes | c.2369C>T | Likely damaging | 0.92 | 0.97 | 0.53 | 27.5 | -3.91 |
| G810D | 3 | 4 | yes | c.2429G>A | Likely damaging | 0.90 | 0.96 | 0.49 | 28.8 | -3.65 |
| L833V | 38 | 2 | yes | c.2497T>G | Uncertain (conflicting predictors) | 0.41 | 0.63 | 0.37 | 18.5 | -3.80 |
| V834L | 12 | 2 | yes | c.2500G>T, c.2500G>C | Likely damaging | 0.99 | 0.99 | 0.82 | 30.4 | -3.88 |
| H835L | 16 | 1 | yes | c.2504A>T | Likely damaging | 1.00 | 1.00 | 0.92 | 31.0 | -4.52 |
| P848L | 26 | 1 | yes | c.2543C>T | Likely damaging | 0.82 | 0.74 | 0.47 | 28.8 | -3.77 |
| T854A | 4 | 1 | yes | c.2560A>G | Likely damaging | 0.92 | 0.76 | 0.70 | 28.1 | -3.53 |
| G857E | 5 | 1 | no | c.2570G>A | Likely damaging | 1.00 | 1.00 | 0.96 | 33.0 | -4.40 |
| L858A | 2 | 3 | yes | c.2572_2573delinsGC | Likely damaging (single predictor) |  |  |  |  |  |
| L858M | 7 | 3 | yes | c.2572C>A | Likely damaging | 0.97 | 0.93 | 0.80 | 24.7 | -3.86 |
| L858R | 10388 | 4 | yes | c.2573T>G | Likely damaging | 1.00 | 1.00 | 0.96 | 33.0 | -4.33 |
| L861P | 2 | 2 | yes | c.2582T>C | Likely damaging | 0.99 | 1.00 | 0.92 | 30.0 | -4.04 |
| L861Q | 418 | 3 | yes | c.2582T>A | Likely damaging | 0.99 | 1.00 | 0.86 | 29.7 | -4.02 |
| L861R | 19 | 2 | yes | c.2582T>G | Likely damaging | 0.99 | 0.99 | 0.89 | 29.8 | -4.00 |
| A864T | 5 | 1 | no | c.2590G>A | Uncertain (conflicting predictors) | 0.41 | 0.14 | 0.45 | 33.0 | -3.16 |
| H870R | 7 | 2 | yes | c.2609A>G | Likely benign | 0.29 | 0.22 | 0.43 | 24.1 | -3.50 |

## In-frame deletions and insertions

Typed in protein notation; the pipeline maps a protein deletion to one DNA change and cannot map a protein
insertion (many DNA sequences give it). In-frame changes get no model score: they are called *likely
damaging* when CADD is at least 20.

87 variants; 21 annotated, 66 not mapped (66 × an insertion written as a protein change cannot be traced to one DNA change: give the c. change). Ensembl's CADD scored 0 of the 21 annotated deletions, so they stay *protein-altering in-frame* (uncertain).

| 3T3 score | Damaging | Uncertain | Benign / silent | Not assessed | mixed | error |
|---|---:|---:|---:|---:|---:|---:|
| 4 | 0 | 4 | 0 | 0 | 0 | 38 |
| 3 | 0 | 4 | 0 | 0 | 0 | 22 |
| 2 | 0 | 3 | 0 | 0 | 0 | 3 |
| 1 | 0 | 10 | 0 | 0 | 0 | 3 |

9 variants with several changes on one allele (such as L858R_T790M) are not single inputs
and are left out.

## Every missense variant

<details><summary>Table</summary>

| Variant | COSMIC | 3T3 | Ba/F3 | DNA change | Call | Model | AlphaMissense | REVEL | CADD | popEVE |
|---|---:|---:|---|---|---|---:|---:|---:|---:|---:|
| R23W | 1 | 1 | yes | c.67C>T | Likely benign | 0.03 | 0.13 | 0.22 | 14.8 | -3.18 |
| H47Y | 6 | 1 | no | c.139C>T | Likely damaging | 0.95 | 0.90 | 0.80 | 24.2 | -3.53 |
| L62R | 24 | 2 | yes | c.185T>G | Likely damaging | 0.93 | 0.74 | 0.80 | 26.2 | -3.46 |
| G63R | 5 | 1 | yes | c.187G>A, c.187G>C | Likely damaging | 0.98 | 0.85 | 0.85 | 29.6 | -3.81 |
| E84V | 3 | 1 | yes | c.251A>T | Likely damaging | 0.98 | 0.98 | 0.85 | 28.5 | -3.83 |
| R98Q | 5 | 2 | yes | c.293G>A | Likely benign | 0.01 | 0.07 | 0.08 | 14.4 | -3.08 |
| R108K | 45 | 3 | yes | c.323G>A | Likely damaging | 0.98 | 0.91 | 0.85 | 28.6 | -3.75 |
| E114K | 13 | 1 | yes | c.340G>A | Uncertain (intermediate scores) | 0.45 | 0.56 | 0.46 | 22.3 | -3.24 |
| E160K | 2 | 1 | yes | c.478G>A | Likely benign | 0.25 | 0.12 | 0.50 | 23.4 | -3.27 |
| R165Q | 4 | 3 | yes | c.494G>A | Likely benign | 0.01 | 0.06 | 0.18 | 13.5 | -2.90 |
| C187R | 2 | 2 | yes | c.559T>C | Likely damaging | 0.85 | 0.88 | 0.61 | 23.0 | -3.41 |
| G203R | 2 | 1 | yes | c.607G>C, c.607G>A | Uncertain (conflicting predictors) | 0.71 | 0.29 | 0.75 | 25.1 | -3.42 |
| R222C | 21 | 3 | yes | c.664C>T | Likely damaging | 0.98 | 0.99 | 0.71 | 32.0 | -4.03 |
| S246R | 2 | 1 | yes | c.738C>A, c.736A>C, c.738C>G | Likely benign / Uncertain (conflicting predictors) | 0.36 | 0.75 | 0.35 | 17.8 | -3.27 |
| R252C | 12 | 1 | yes | c.754C>T | Likely damaging | 0.92 | 0.72 | 0.58 | 32.0 | -3.88 |
| R252H | 8 | 2 | yes | c.755G>A | Likely benign | 0.19 | 0.15 | 0.40 | 23.1 | -3.37 |
| R252P | 3 | 1 | yes | c.755G>C | Likely damaging | 0.96 | 0.96 | 0.78 | 24.9 | -3.68 |
| D256G | 3 | 2 | yes | c.767A>G | Likely damaging | 0.93 | 0.93 | 0.63 | 27.8 | -3.62 |
| T263P | 18 | 2 | yes | c.787A>C | Uncertain (conflicting predictors) | 0.66 | 0.47 | 0.55 | 27.4 | -3.33 |
| Y275F | 2 | 2 | yes | c.824A>T | Likely benign | 0.08 | 0.22 | 0.21 | 23.7 | -2.94 |
| A289D | 4 | 1 | yes | c.866C>A | Likely damaging | 0.94 | 0.95 | 0.59 | 29.4 | -3.70 |
| A289T | 31 | 1 | yes | c.865G>A | Likely damaging | 0.85 | 0.61 | 0.57 | 31.0 | -3.61 |
| A289V | 12 | 1 | yes | c.866C>T | Likely damaging | 0.94 | 0.83 | 0.65 | 29.4 | -3.75 |
| H304Y | 6 | 1 | yes | c.910C>T | Likely benign | 0.13 | 0.21 | 0.28 | 23.6 | -3.27 |
| G312W | 3 | 1 | yes | c.934G>T | Likely benign | 0.19 | 0.34 | 0.35 | 20.1 | -3.37 |
| V323I | 4 | 1 | yes | c.967G>A | Likely benign | 0.02 | 0.06 | 0.09 | 18.2 | -2.97 |
| R324L | 8 | 1 | yes | c.971G>T | Uncertain (conflicting predictors) | 0.70 | 0.29 | 0.69 | 25.8 | -3.54 |
| G331R | 2 | 1 | yes | c.991G>A, c.991G>C | Likely damaging | 0.82 | 0.81 | 0.45 | 29.2 | -3.60 |
| T363I | 8 | 2 | yes | c.1088C>T | Likely damaging | 0.96 | 0.71 | 0.89 | 26.4 | -3.80 |
| T384S | 2 | 1 | yes | c.1150A>T, c.1151C>G | Likely benign | 0.18 | 0.45 | 0.23 | 23.0 | -3.27 |
| R427C | 4 | 1 | yes | c.1279C>T | Likely damaging | 0.99 | 0.96 | 0.80 | 33.0 | -3.96 |
| I462M | 4 | 1 | yes | c.1386A>G | Likely benign | 0.04 | 0.19 | 0.40 | 7.9 | -3.06 |
| S492R | 8 | 1 | yes | c.1476C>A, c.1474A>C, c.1476C>G | Likely benign | 0.03 | 0.19 | 0.17 | 15.4 | -2.96 |
| R521K | 120 | 2 | yes | c.1562G>A | Likely benign | 0.00 | 0.06 | 0.01 | 0.2 | -2.89 |
| G573_STOP | 3 | 1 | no | c.1717G>T | Loss of function (truncating) |  |  |  | 40.0 |  |
| P589L | 4 | 2 | yes | c.1766C>T | Likely benign | 0.08 | 0.23 | 0.18 | 22.1 | -3.40 |
| P596L | 16 | 2 | yes | c.1787C>T | Likely damaging | 0.92 | 0.76 | 0.64 | 29.0 | -3.80 |
| P596R | 2 | 2 | yes | c.1787C>G | Likely damaging | 0.94 | 0.90 | 0.64 | 27.0 | -3.85 |
| P596S | 2 | 1 | yes | c.1786C>T | Likely damaging | 0.81 | 0.80 | 0.49 | 26.0 | -3.70 |
| G598A | 5 | 1 | yes | c.1793G>C | Uncertain (conflicting predictors) | 0.42 | 0.36 | 0.44 | 25.4 | -3.50 |
| G598V | 69 | 4 | yes | c.1793G>T | Likely damaging | 0.90 | 0.84 | 0.61 | 26.5 | -3.78 |
| N604D | 2 | 2 | yes | c.1810A>G | Likely benign | 0.01 | 0.07 | 0.09 | 14.3 | -2.74 |
| T605A | 2 | 1 | yes | c.1813A>G | Likely benign | 0.04 | 0.11 | 0.16 | 21.1 | -3.24 |
| C620Y | 4 | 1 | yes | c.1859G>A | Likely damaging | 0.98 | 0.99 | 0.73 | 29.7 | -4.11 |
| G627E | 2 | 1 | yes | c.1880G>A | Uncertain (intermediate scores) | 0.66 | 0.53 | 0.39 | 31.0 | -3.67 |
| C628F | 2 | 1 | yes | c.1883G>T | Likely damaging | 0.93 | 0.97 | 0.59 | 24.6 | -4.02 |
| C636Y | 3 | 1 | yes | c.1907G>A | Likely damaging | 0.91 | 0.95 | 0.55 | 26.1 | -3.83 |
| T638M | 4 | 4 | yes | c.1913C>T | Likely benign | 0.01 | 0.08 | 0.27 | 0.0 | -3.15 |
| P641H | 2 | 4 | yes | c.1922C>A | Likely benign | 0.04 | 0.10 | 0.30 | 14.7 | -3.13 |
| I643V | 2 | 1 | yes | c.1927A>G | Likely benign | 0.01 | 0.12 | 0.11 | 13.7 | -2.86 |
| S645C | 6 | 3 | yes | c.1934C>G | Uncertain (intermediate scores) | 0.44 | 0.38 | 0.47 | 23.5 | -3.62 |
| V651M | 4 | 3 | yes | c.1951G>A | Uncertain (intermediate scores) | 0.49 | 0.27 | 0.57 | 24.0 | -3.54 |
| R671C | 3 | 2 | yes | c.2011C>T | Likely damaging | 0.83 | 0.21 | 0.73 | 32.0 | -3.75 |
| R677H | 6 | 3 | yes | c.2030G>A | Uncertain (conflicting predictors) | 0.46 | 0.41 | 0.44 | 25.5 | -3.47 |
| R680Q | 3 | 3 | yes | c.2039G>A | Likely damaging | 0.89 | 0.55 | 0.66 | 31.0 | -3.66 |
| L688P | 2 | 1 | no | c.2063T>C | Likely damaging | 0.96 | 0.98 | 0.66 | 27.1 | -3.95 |
| V689L | 2 | 1 | yes | c.2065G>T, c.2065G>C | Uncertain (conflicting predictors) | 0.66 | 0.81 | 0.38 | 25.4 | -3.43 |
| P691S | 3 | 1 | yes | c.2071C>T | Likely damaging | 0.94 | 0.98 | 0.61 | 26.6 | -3.82 |
| L692P | 3 | 1 | no | c.2075T>C | Likely damaging | 0.94 | 0.98 | 0.60 | 27.0 | -3.96 |
| P694L | 3 | 1 | yes | c.2081C>T | Likely damaging | 0.93 | 0.94 | 0.60 | 27.3 | -3.76 |
| P694S | 5 | 1 | yes | c.2080C>T | Likely damaging | 0.91 | 0.97 | 0.56 | 25.9 | -3.72 |
| A698T | 2 | 1 | yes | c.2092G>A | Uncertain (intermediate scores) | 0.31 | 0.30 | 0.39 | 24.4 | -3.53 |
| P699S | 2 | 1 | no | c.2095C>T | Likely damaging | 0.93 | 0.96 | 0.60 | 25.3 | -3.84 |
| A702D | 2 | 1 | no | c.2105C>A | Likely damaging | 0.95 | 0.98 | 0.65 | 26.5 | -3.99 |
| A702S | 2 | 1 | yes | c.2104G>T | Uncertain (conflicting predictors) | 0.51 | 0.37 | 0.49 | 25.4 | -3.60 |
| L703P | 4 | 1 | no | c.2108T>C | Likely damaging | 0.87 | 0.89 | 0.61 | 24.2 | -3.52 |
| L703V | 5 | 2 | yes | c.2107C>G | Likely benign | 0.12 | 0.22 | 0.24 | 22.8 | -3.56 |
| R705K | 3 | 1 | yes | c.2114G>A | Likely damaging | 0.89 | 0.78 | 0.62 | 26.8 | -3.65 |
| I706T | 4 | 3 | yes | c.2117T>C | Likely damaging | 0.98 | 0.96 | 0.85 | 26.5 | -3.84 |
| K708E | 2 | 1 | no | c.2122A>G | Likely damaging | 0.93 | 0.99 | 0.51 | 27.7 | -4.12 |
| K708N | 2 | 1 | yes | c.2124G>C, c.2124G>T | Likely damaging | 0.76 | 0.99 | 0.30 | 25.6 | -3.88 |
| E709A | 37 | 3 | yes | c.2126A>C | Likely damaging | 0.93 | 0.90 | 0.59 | 27.2 | -3.91 |
| E709G | 18 | 2 | yes | c.2126A>G | Likely damaging | 0.96 | 0.87 | 0.67 | 29.0 | -4.00 |
| E709H | 2 | 3 | yes | c.2125_2127delinsCAT, c.2125_2127delinsCAC | Likely damaging (single predictor) |  |  |  |  |  |
| E709K | 63 | 4 | yes | c.2125G>A | Likely damaging | 0.92 | 0.96 | 0.47 | 29.1 | -4.04 |
| E709V | 8 | 4 | yes | c.2126A>T | Likely damaging | 0.97 | 0.96 | 0.58 | 34.0 | -4.03 |
| F712L | 2 | 1 | yes | c.2136C>A, c.2134T>C, c.2136C>G | Likely benign | 0.19 | 0.75 | 0.26 | 17.3 | -2.88 |
| F712S | 2 | 1 | yes | c.2135T>C | Likely damaging | 0.98 | 0.99 | 0.77 | 28.5 | -4.04 |
| K713F | 2 | 1 | yes | c.2137_2139delinsTTC, c.2137_2139inv | Uncertain (intermediate scores) |  |  |  |  |  |
| K714N | 6 | 1 | yes | c.2142G>C, c.2142G>T | Likely damaging | 0.86 | 0.97 | 0.52 | 23.4 | -3.79 |
| L718M | 4 | 1 | yes | c.2152C>A | Likely damaging | 0.84 | 0.77 | 0.58 | 22.4 | -4.02 |
| L718P | 3 | 1 | no | c.2153T>C | Likely damaging | 0.99 | 1.00 | 0.89 | 28.0 | -4.41 |
| L718V | 2 | 1 | yes | c.2152C>G | Likely damaging | 0.80 | 0.79 | 0.54 | 21.9 | -3.88 |
| G719A | 206 | 3 | yes | c.2156G>C | Likely damaging | 0.98 | 0.99 | 0.82 | 25.9 | -4.04 |
| G719C | 93 | 4 | yes | c.2155G>T | Likely damaging | 0.99 | 1.00 | 0.88 | 28.9 | -4.30 |
| G719D | 16 | 3 | yes | c.2156G>A | Likely damaging | 0.99 | 1.00 | 0.86 | 26.9 | -4.24 |
| G719S | 159 | 3 | yes | c.2155G>A | Likely damaging | 0.99 | 1.00 | 0.85 | 27.6 | -4.19 |
| S720F | 16 | 2 | yes | c.2159C>T | Likely damaging | 0.85 | 0.66 | 0.60 | 25.5 | -3.97 |
| S720P | 4 | 1 | yes | c.2158T>C | Likely damaging | 0.96 | 0.99 | 0.65 | 26.1 | -4.16 |
| G721A | 2 | 1 | no | c.2162G>C | Likely damaging | 0.98 | 1.00 | 0.84 | 25.4 | -4.27 |
| G721S | 3 | 1 | no | c.2161G>A | Likely damaging | 0.99 | 1.00 | 0.86 | 27.6 | -4.23 |
| A722V | 2 | 1 | yes | c.2165C>T | Likely damaging | 0.87 | 0.98 | 0.49 | 24.3 | -3.87 |
| G724S | 15 | 1 | yes | c.2170G>A | Likely damaging | 0.97 | 1.00 | 0.74 | 27.6 | -3.98 |
| T725M | 10 | 1 | yes | c.2174C>T | Likely damaging | 0.80 | 0.40 | 0.68 | 26.5 | -3.83 |
| Y727C | 3 | 1 | yes | c.2180A>G | Likely damaging | 0.86 | 0.43 | 0.74 | 26.9 | -3.86 |
| Y727H | 2 | 2 | yes | c.2179T>C | Uncertain (intermediate scores) | 0.40 | 0.32 | 0.47 | 25.2 | -3.39 |
| G729E | 3 | 1 | no | c.2186G>A | Likely damaging | 1.00 | 1.00 | 0.96 | 34.0 | -4.25 |
| G729R | 2 | 1 | no | c.2185G>C, c.2185G>A | Likely damaging | 1.00 | 1.00 | 0.95 | 35.5 | -4.29 |
| L730F | 3 | 1 | yes | c.2188C>T | Likely benign | 0.14 | 0.31 | 0.20 | 24.6 | -3.42 |
| W731_STOP | 7 | 1 | no | c.2193G>A, c.2192G>A | Loss of function (truncating) |  |  |  | 51.5 |  |
| P733L | 4 | 1 | no | c.2198C>T | Likely damaging | 0.98 | 0.97 | 0.78 | 32.0 | -3.87 |
| P733S | 7 | 1 | yes | c.2197C>T | Likely damaging | 0.98 | 1.00 | 0.74 | 29.0 | -3.86 |
| E734K | 3 | 3 | yes | c.2200G>A | Likely damaging | 0.93 | 0.83 | 0.52 | 34.0 | -3.86 |
| G735S | 7 | 1 | yes | c.2203G>A | Likely damaging | 0.86 | 0.55 | 0.56 | 32.0 | -3.82 |
| V738G | 2 | 1 | no | c.2213T>G | Likely damaging | 0.98 | 0.92 | 0.77 | 29.7 | -3.97 |
| P741L | 5 | 1 | no | c.2222C>T | Likely damaging | 0.95 | 0.74 | 0.73 | 29.1 | -3.93 |
| P741S | 3 | 1 | yes | c.2221C>T | Likely damaging | 0.90 | 0.60 | 0.75 | 27.8 | -3.55 |
| V742A | 6 | 1 | yes | c.2225T>C | Likely damaging | 0.99 | 0.97 | 0.90 | 27.2 | -3.86 |
| V742I | 7 | 1 | yes | c.2224G>A | Likely damaging | 0.84 | 0.39 | 0.65 | 31.0 | -3.81 |
| A743S | 2 | 1 | yes | c.2227G>T | Likely damaging | 0.98 | 0.84 | 0.81 | 32.0 | -4.08 |
| A743T | 4 | 1 | no | c.2227G>A | Likely damaging | 0.99 | 0.99 | 0.91 | 33.0 | -4.25 |
| A743V | 4 | 1 | no | c.2228C>T | Likely damaging | 0.99 | 0.99 | 0.91 | 29.2 | -4.13 |
| E746K | 3 | 1 | yes | c.2236G>A | Uncertain (conflicting predictors) | 0.68 | 0.88 | 0.27 | 27.2 | -3.80 |
| E746V | 5 | 1 | yes | c.2237A>T | Likely benign | 0.05 | 0.27 | 0.20 | 21.4 | -2.70 |
| L747P | 27 | 4 | yes | c.2239_2240delinsCC | Likely damaging (single predictor) |  |  |  |  |  |
| L747S | 16 | 4 | yes | c.2240T>C | Likely damaging | 0.99 | 1.00 | 0.84 | 27.3 | -4.20 |
| R748K | 4 | 1 | yes | c.2243G>A | Uncertain (intermediate scores) | 0.54 | 0.52 | 0.43 | 26.7 | -3.37 |
| E749G | 2 | 1 | yes | c.2246A>G | Likely damaging | 0.96 | 0.87 | 0.75 | 28.0 | -3.88 |
| A750P | 9 | 1 | yes | c.2248G>C | Likely damaging | 0.91 | 0.87 | 0.62 | 27.5 | -3.60 |
| T751I | 10 | 1 | no | c.2252C>T | Likely damaging | 0.91 | 0.81 | 0.65 | 26.5 | -3.77 |
| S752P | 2 | 1 | yes | c.2254T>C | Likely damaging | 0.97 | 0.94 | 0.79 | 26.3 | -3.89 |
| P753L | 9 | 1 | yes | c.2258C>T | Likely damaging | 0.76 | 0.53 | 0.58 | 27.3 | -3.62 |
| P753S | 18 | 1 | yes | c.2257C>T | Uncertain (intermediate scores) | 0.49 | 0.39 | 0.51 | 24.9 | -3.35 |
| K754E | 6 | 1 | yes | c.2260A>G | Likely damaging | 0.89 | 0.74 | 0.64 | 26.6 | -3.81 |
| K754Q | 2 | 1 | yes | c.2260A>C | Uncertain (conflicting predictors) | 0.46 | 0.47 | 0.39 | 25.8 | -3.59 |
| A755D | 2 | 1 | yes | c.2264C>A | Likely damaging | 0.98 | 0.96 | 0.81 | 28.2 | -3.76 |
| K757R | 3 | 1 | yes | c.2270A>G | Likely benign | 0.23 | 0.09 | 0.46 | 25.3 | -3.19 |
| D761G | 2 | 1 | yes | c.2282A>G | Likely damaging | 0.99 | 0.90 | 0.86 | 36.0 | -3.88 |
| D761N | 9 | 1 | yes | c.2281G>A | Likely damaging | 0.91 | 0.72 | 0.53 | 34.0 | -3.87 |
| D761Y | 12 | 1 | yes | c.2281G>T | Likely damaging | 0.99 | 0.95 | 0.77 | 34.0 | -4.22 |
| V765M | 11 | 1 | yes | c.2293G>A | Uncertain (conflicting predictors) | 0.65 | 0.68 | 0.29 | 32.0 | -3.50 |
| A767V | 6 | 2 | yes | c.2300C>T | Likely damaging | 0.99 | 0.99 | 0.78 | 32.0 | -3.98 |
| S768I | 273 | 4 | yes | c.2303G>T | Likely damaging | 0.98 | 0.90 | 0.77 | 31.0 | -3.96 |
| V769L | 17 | 2 | yes | c.2305G>C, c.2305G>T | Uncertain (conflicting predictors) | 0.59 | 0.86 | 0.29 | 25.6 | -3.44 |
| V769M | 12 | 2 | yes | c.2305G>A | Uncertain (conflicting predictors) | 0.61 | 0.57 | 0.32 | 33.0 | -3.40 |
| D770N | 5 | 4 | yes | c.2308G>A | Likely benign | 0.07 | 0.16 | 0.12 | 25.8 | -3.30 |
| P772S | 2 | 1 | yes | c.2314C>T | Uncertain (conflicting predictors) | 0.68 | 0.26 | 0.63 | 28.2 | -3.66 |
| H773L | 11 | 2 | yes | c.2318A>T | Likely damaging | 0.88 | 0.53 | 0.66 | 28.4 | -4.08 |
| H773P | 4 | 1 | no | c.2318A>C | Likely damaging | 0.99 | 0.92 | 0.77 | 32.0 | -4.31 |
| H773R | 11 | 1 | yes | c.2318A>G | Likely damaging | 0.90 | 0.83 | 0.59 | 25.1 | -4.06 |
| H773Y | 5 | 1 | yes | c.2317C>T | Likely benign | 0.08 | 0.13 | 0.23 | 24.1 | -3.12 |
| V774M | 23 | 3 | yes | c.2320G>A | Likely damaging | 0.97 | 0.93 | 0.62 | 33.0 | -3.96 |
| C775Y | 4 | 1 | no | c.2324G>A | Likely damaging | 0.97 | 0.99 | 0.58 | 33.0 | -4.07 |
| R776C | 9 | 1 | yes | c.2326C>T | Likely damaging | 0.96 | 0.84 | 0.59 | 34.0 | -3.92 |
| R776H | 21 | 1 | yes | c.2327G>A | Likely damaging | 0.95 | 0.65 | 0.68 | 34.0 | -3.81 |
| L777Q | 2 | 1 | yes | c.2330T>A | Likely damaging | 0.99 | 1.00 | 0.86 | 33.0 | -4.33 |
| G779C | 3 | 4 | yes | c.2335G>T | Likely damaging | 1.00 | 1.00 | 0.94 | 33.0 | -4.19 |
| G779D | 3 | 1 | no | c.2336G>A | Likely damaging | 1.00 | 1.00 | 0.92 | 33.0 | -4.25 |
| G779F | 4 | 1 | yes | c.2335_2336delinsTT | Likely damaging (single predictor) |  |  |  |  |  |
| G779S | 6 | 1 | yes | c.2335G>A | Likely damaging | 1.00 | 0.99 | 0.94 | 33.0 | -4.24 |
| T783I | 4 | 1 | yes | c.2348C>T | Likely damaging | 0.97 | 0.98 | 0.76 | 27.2 | -3.98 |
| S784F | 7 | 1 | yes | c.2351C>T | Likely damaging | 0.90 | 0.95 | 0.46 | 27.7 | -4.07 |
| S784P | 3 | 1 | yes | c.2350T>C | Likely benign | 0.27 | 0.64 | 0.21 | 23.7 | -3.25 |
| T785I | 6 | 1 | yes | c.2354C>T | Uncertain (intermediate scores) | 0.67 | 0.66 | 0.42 | 25.8 | -3.76 |
| V786M | 12 | 1 | yes | c.2356G>A | Uncertain (intermediate scores) | 0.70 | 0.79 | 0.37 | 27.3 | -3.61 |
| L788F | 2 | 1 | yes | c.2362C>T | Likely damaging | 0.99 | 0.99 | 0.86 | 27.5 | -4.04 |
| T790A | 2 | 1 | yes | c.2368A>G | Likely damaging | 0.96 | 0.96 | 0.69 | 26.3 | -3.97 |
| T790M | 1398 | 3 | yes | c.2369C>T | Likely damaging | 0.92 | 0.97 | 0.53 | 27.5 | -3.91 |
| L792P | 2 | 1 | no | c.2375T>C | Likely damaging | 0.98 | 1.00 | 0.74 | 28.8 | -4.18 |
| F795S | 2 | 1 | no | c.2384T>C | Uncertain (intermediate scores) | 0.42 | 0.63 | 0.33 | 23.5 | -3.48 |
| G796A | 3 | 1 | no | c.2387G>C | Likely damaging | 0.99 | 0.98 | 0.95 | 28.4 | -4.23 |
| G796D | 6 | 1 | no | c.2387G>A | Likely damaging | 0.99 | 1.00 | 0.94 | 29.2 | -4.32 |
| G796S | 8 | 1 | no | c.2386G>A | Likely damaging | 0.99 | 0.95 | 0.92 | 33.0 | -4.18 |
| C797S | 11 | 1 | no | c.2390G>C, c.2389T>A | Likely damaging / Uncertain (conflicting predictors) | 0.77 | 0.74 | 0.48 | 27.2 | -3.58 |
| C797Y | 3 | 1 | no | c.2390G>A | Likely damaging | 0.97 | 0.99 | 0.62 | 32.0 | -4.20 |
| L798H | 3 | 1 | no | c.2393T>A | Likely damaging | 0.99 | 1.00 | 0.94 | 28.3 | -4.46 |
| L799R | 2 | 4 | yes | c.2396T>G | Likely damaging | 0.94 | 0.95 | 0.61 | 26.1 | -4.03 |
| Y801C | 2 | 1 | no | c.2402A>G | Likely damaging | 0.99 | 0.90 | 0.94 | 26.5 | -4.03 |
| Y801H | 2 | 4 | yes | c.2401T>C | Likely damaging | 0.94 | 0.75 | 0.78 | 26.2 | -3.63 |
| V802A | 2 | 1 | no | c.2405T>C | Likely damaging | 0.95 | 0.90 | 0.74 | 25.7 | -3.81 |
| V802F | 3 | 1 | yes | c.2404G>T | Likely damaging | 0.95 | 0.95 | 0.70 | 23.5 | -4.15 |
| V802I | 5 | 1 | yes | c.2404G>A | Likely benign | 0.12 | 0.09 | 0.49 | 17.0 | -3.25 |
| R803L | 2 | 1 | no | c.2408G>T | Likely damaging | 0.88 | 0.67 | 0.65 | 26.5 | -3.83 |
| E804G | 2 | 1 | yes | c.2411A>G | Likely damaging | 0.81 | 0.31 | 0.75 | 27.3 | -3.77 |
| E804K | 3 | 1 | yes | c.2410G>A | Uncertain (intermediate scores) | 0.55 | 0.24 | 0.62 | 25.1 | -3.53 |
| H805L | 2 | 2 | yes | c.2414A>T | Likely damaging | 0.95 | 0.68 | 0.81 | 26.5 | -3.92 |
| K806E | 4 | 1 | yes | c.2416A>G | Uncertain (conflicting predictors) | 0.54 | 0.33 | 0.44 | 27.9 | -3.85 |
| K806R | 3 | 1 | yes | c.2417A>G | Likely benign | 0.05 | 0.08 | 0.17 | 21.7 | -3.26 |
| G810A | 2 | 1 | yes | c.2429G>C | Likely damaging | 0.76 | 0.66 | 0.47 | 28.0 | -3.71 |
| G810D | 3 | 4 | yes | c.2429G>A | Likely damaging | 0.90 | 0.96 | 0.49 | 28.8 | -3.65 |
| G810S | 3 | 1 | yes | c.2428G>A | Uncertain (conflicting predictors) | 0.52 | 0.36 | 0.41 | 29.5 | -3.61 |
| S811F | 3 | 2 | yes | c.2432C>T | Likely damaging | 0.95 | 1.00 | 0.56 | 29.5 | -3.99 |
| L814P | 3 | 1 | yes | c.2441T>C | Likely damaging | 0.99 | 1.00 | 0.79 | 28.8 | -4.27 |
| V819A | 2 | 1 | yes | c.2456T>C | Likely damaging | 0.92 | 0.92 | 0.56 | 28.0 | -3.78 |
| I821T | 2 | 1 | no | c.2462T>C | Likely damaging | 0.97 | 0.96 | 0.77 | 26.6 | -4.01 |
| A822P | 1 | 1 | no | c.2464G>C | Likely damaging | 0.99 | 1.00 | 0.91 | 29.5 | -4.45 |
| N826D | 2 | 1 | no | c.2476A>G | Uncertain (intermediate scores) | 0.40 | 0.41 | 0.41 | 24.7 | -3.53 |
| E829K | 2 | 1 | yes | c.2485G>A | Likely damaging | 0.98 | 0.99 | 0.66 | 33.0 | -4.22 |
| R831C | 5 | 1 | no | c.2491C>T | Likely damaging | 0.86 | 0.19 | 0.77 | 33.0 | -3.67 |
| R831H | 10 | 1 | yes | c.2492G>A | Uncertain (conflicting predictors) | 0.42 | 0.12 | 0.56 | 27.1 | -3.36 |
| R831L | 3 | 1 | yes | c.2492G>T | Likely damaging | 0.85 | 0.54 | 0.66 | 29.3 | -3.57 |
| R832C | 4 | 1 | yes | c.2494C>T | Uncertain (conflicting predictors) | 0.70 | 0.37 | 0.49 | 33.0 | -3.62 |
| R832H | 4 | 1 | yes | c.2495G>A | Likely benign | 0.09 | 0.13 | 0.23 | 25.1 | -3.19 |
| L833F | 7 | 3 | yes | c.2499G>T, c.2499G>C | Likely damaging | 0.80 | 0.97 | 0.37 | 25.9 | -3.74 |
| L833V | 38 | 2 | yes | c.2497T>G | Uncertain (conflicting predictors) | 0.41 | 0.63 | 0.37 | 18.5 | -3.80 |
| V834A | 3 | 1 | no | c.2501T>C | Likely damaging | 0.99 | 0.98 | 0.87 | 29.7 | -4.04 |
| V834L | 12 | 2 | yes | c.2500G>T, c.2500G>C | Likely damaging | 0.99 | 0.99 | 0.82 | 30.4 | -3.88 |
| V834M | 3 | 1 | yes | c.2500G>A | Likely damaging | 0.99 | 0.99 | 0.84 | 32.0 | -4.09 |
| H835L | 16 | 1 | yes | c.2504A>T | Likely damaging | 1.00 | 1.00 | 0.92 | 31.0 | -4.52 |
| R836C | 5 | 1 | yes | c.2506C>T | Likely damaging | 0.99 | 0.99 | 0.88 | 33.0 | -4.27 |
| R836H | 2 | 1 | no | c.2507G>A | Likely damaging | 0.99 | 1.00 | 0.86 | 34.0 | -4.36 |
| R836S | 3 | 1 | no | c.2506C>A | Likely damaging | 0.99 | 1.00 | 0.87 | 29.5 | -4.39 |
| D837N | 8 | 1 | no | c.2509G>A | Likely damaging | 0.98 | 1.00 | 0.67 | 33.0 | -3.93 |
| L838M | 2 | 2 | yes | c.2512C>A | Likely damaging | 0.94 | 0.95 | 0.66 | 24.7 | -4.05 |
| L838P | 4 | 1 | no | c.2513T>C | Likely damaging | 1.00 | 1.00 | 0.91 | 33.0 | -4.46 |
| A839T | 6 | 1 | no | c.2515G>A | Likely damaging | 0.98 | 1.00 | 0.65 | 33.0 | -4.22 |
| A840T | 5 | 1 | yes | c.2518G>A | Likely damaging | 0.96 | 0.96 | 0.60 | 32.0 | -3.86 |
| R841K | 2 | 1 | no | c.2522G>A | Likely damaging | 0.99 | 0.99 | 0.82 | 33.0 | -4.37 |
| N842D | 2 | 1 | no | c.2524A>G | Likely damaging | 1.00 | 1.00 | 0.90 | 32.0 | -4.53 |
| V843I | 15 | 3 | yes | c.2527G>A | Likely benign | 0.23 | 0.18 | 0.25 | 31.0 | -3.46 |
| T847A | 3 | 1 | yes | c.2539A>G | Uncertain (conflicting predictors) | 0.46 | 0.26 | 0.51 | 26.1 | -3.49 |
| T847I | 8 | 1 | no | c.2540C>T | Likely damaging | 0.91 | 0.85 | 0.59 | 26.9 | -3.90 |
| P848L | 26 | 1 | yes | c.2543C>T | Likely damaging | 0.82 | 0.74 | 0.47 | 28.8 | -3.77 |
| Q849R | 2 | 1 | yes | c.2546A>G | Likely benign | 0.22 | 0.30 | 0.28 | 26.0 | -3.42 |
| H850D | 5 | 1 | no | c.2548C>G | Likely damaging | 0.94 | 0.94 | 0.54 | 29.5 | -4.06 |
| V851I | 7 | 1 | yes | c.2551G>A | Likely benign | 0.25 | 0.12 | 0.49 | 24.0 | -3.30 |
| T854A | 4 | 1 | yes | c.2560A>G | Likely damaging | 0.92 | 0.76 | 0.70 | 28.1 | -3.53 |
| D855N | 4 | 1 | no | c.2563G>A | Likely damaging | 0.99 | 1.00 | 0.85 | 33.0 | -4.21 |
| F856L | 3 | 1 | no | c.2568T>A, c.2566T>C, c.2568T>G | Likely damaging | 0.96 | 1.00 | 0.82 | 23.5 | -3.97 |
| G857E | 5 | 1 | no | c.2570G>A | Likely damaging | 1.00 | 1.00 | 0.96 | 33.0 | -4.40 |
| G857R | 5 | 1 | yes | c.2569G>A, c.2569G>C | Likely damaging | 1.00 | 1.00 | 0.95 | 33.5 | -4.32 |
| L858A | 2 | 3 | yes | c.2572_2573delinsGC | Likely damaging (single predictor) |  |  |  |  |  |
| L858M | 7 | 3 | yes | c.2572C>A | Likely damaging | 0.97 | 0.93 | 0.80 | 24.7 | -3.86 |
| L858P | 2 | 1 | no | c.2573T>C | Likely damaging | 1.00 | 1.00 | 0.96 | 33.0 | -4.31 |
| L858R | 10388 | 4 | yes | c.2573T>G | Likely damaging | 1.00 | 1.00 | 0.96 | 33.0 | -4.33 |
| L858W | 2 | 1 | no | c.2572_2573delinsTG | Likely damaging (single predictor) |  |  |  |  |  |
| A859T | 12 | 1 | yes | c.2575G>A | Likely damaging | 0.98 | 0.96 | 0.75 | 32.0 | -3.70 |
| K860I | 15 | 2 | yes | c.2579A>T | Likely damaging | 0.99 | 0.99 | 0.79 | 31.0 | -4.10 |
| L861F | 2 | 3 | yes | c.2581_2583delinsTTT, c.2581_2583delinsTTC | Likely damaging (single predictor) |  |  |  |  |  |
| L861P | 2 | 2 | yes | c.2582T>C | Likely damaging | 0.99 | 1.00 | 0.92 | 30.0 | -4.04 |
| L861Q | 418 | 3 | yes | c.2582T>A | Likely damaging | 0.99 | 1.00 | 0.86 | 29.7 | -4.02 |
| L861R | 19 | 2 | yes | c.2582T>G | Likely damaging | 0.99 | 0.99 | 0.89 | 29.8 | -4.00 |
| G863D | 8 | 2 | yes | c.2588G>A | Uncertain (conflicting predictors) | 0.37 | 0.61 | 0.43 | 22.5 | -2.76 |
| G863S | 4 | 1 | yes | c.2587G>A | Likely benign | 0.16 | 0.19 | 0.40 | 22.7 | -3.02 |
| A864E | 2 | 1 | yes | c.2591C>A | Likely damaging | 0.83 | 0.58 | 0.53 | 34.0 | -3.44 |
| A864T | 5 | 1 | no | c.2590G>A | Uncertain (conflicting predictors) | 0.41 | 0.14 | 0.45 | 33.0 | -3.16 |
| A864V | 6 | 1 | yes | c.2591C>T | Likely benign | 0.26 | 0.20 | 0.46 | 25.2 | -3.11 |
| E866K | 5 | 1 | yes | c.2596G>A | Likely damaging | 0.97 | 0.89 | 0.68 | 33.0 | -3.71 |
| E868G | 4 | 1 | yes | c.2603A>G | Likely damaging | 0.90 | 0.61 | 0.66 | 29.0 | -3.91 |
| E868K | 3 | 1 | no | c.2602G>A | Likely damaging | 0.84 | 0.72 | 0.50 | 28.7 | -3.80 |
| H870R | 7 | 2 | yes | c.2609A>G | Likely benign | 0.29 | 0.22 | 0.43 | 24.1 | -3.50 |
| H870Y | 2 | 1 | yes | c.2608C>T | Uncertain (intermediate scores) | 0.42 | 0.35 | 0.44 | 24.8 | -3.62 |
| A871G | 17 | 3 | yes | c.2612C>G | Likely damaging | 0.90 | 0.73 | 0.59 | 29.3 | -3.82 |
| A871T | 3 | 1 | no | c.2611G>A | Likely damaging | 0.87 | 0.68 | 0.52 | 32.0 | -3.72 |
| A871V | 3 | 1 | yes | c.2612C>T | Likely damaging | 0.93 | 0.82 | 0.58 | 32.0 | -3.72 |
| G873E | 6 | 1 | yes | c.2618G>A | Likely damaging | 0.98 | 0.96 | 0.72 | 32.0 | -3.73 |
| G874S | 5 | 1 | no | c.2620G>A | Likely damaging | 0.93 | 0.87 | 0.52 | 32.0 | -3.85 |
| E884K | 4 | 1 | no | c.2650G>A | Likely damaging | 0.99 | 1.00 | 0.86 | 32.0 | -4.29 |
| G930E | 2 | 1 | no | c.2789G>A | Likely damaging | 0.99 | 1.00 | 0.78 | 31.0 | -4.04 |
| E967K | 4 | 1 | no | c.2899G>A | Likely damaging | 0.96 | 0.99 | 0.53 | 34.0 | -3.96 |
| G983R | 2 | 1 | yes | c.2947G>A, c.2947G>C | Likely damaging | 0.81 | 0.88 | 0.32 | 31.5 | -3.73 |
| D1009N | 7 | 1 | yes | c.3025G>A | Likely benign | 0.18 | 0.13 | 0.35 | 25.3 | -3.41 |
| R1031Q | 3 | 1 | yes | c.3092G>A | Likely benign | 0.12 | 0.09 | 0.30 | 25.9 | -3.19 |
| S1045F | 3 | 1 | yes | c.3134C>T | Likely benign | 0.14 | 0.17 | 0.32 | 23.0 | -3.39 |
| A1048T | 1 | 2 | yes | c.3142G>A | Likely benign | 0.00 | 0.06 | 0.10 | 6.4 | -2.86 |
| A1048V | 2 | 1 | yes | c.3143C>T | Likely benign | 0.01 | 0.09 | 0.09 | 15.2 | -3.00 |
| V1097I | 2 | 1 | yes | c.3289G>A | Likely benign | 0.00 | 0.06 | 0.15 | 1.9 | -2.95 |
| A1118T | 3 | 1 | yes | c.3352G>A | Likely benign | 0.00 | 0.06 | 0.10 | 5.6 | -2.85 |
| S1153I | 3 | 1 | yes | c.3458G>T | Likely benign | 0.01 | 0.08 | 0.08 | 9.5 | -3.26 |
| D1168N | 2 | 1 | yes | c.3502G>A | Uncertain (intermediate scores) | 0.33 | 0.35 | 0.35 | 27.1 | -3.36 |

</details>
