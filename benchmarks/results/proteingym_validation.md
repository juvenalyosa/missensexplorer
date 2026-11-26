# Validation on deep mutational scanning data (ProteinGym)

28701 single-nucleotide missense variants from 71 human assays on 61 proteins of [ProteinGym](https://marks.hms.harvard.edu/proteingym/ProteinGym_v1.3/DMS_ProteinGym_substitutions.zip) (up to 500 per assay). No ClinVar label is involved: each variant's effect was measured in the laboratory. Produced by `benchmarks/validate_proteingym.py`.

- **Spearman**: correlation between the tool's score and the measured loss of function, averaged over the assays of each protein and then over proteins. 95% bootstrap interval over proteins.
- **Δ vs ours**: the tool's mean minus MissensExplorer's, both computed within each assay on the variants both score (paired by variant); an interval that excludes 0 is a real difference.
- **AUC**: against ProteinGym's binary label (fit/unfit), same averaging.
- DMS assays measure one molecular function in one experimental system, not disease, so every tool's correlation is moderate. VARITY used some DMS data in training.

| Tool | Clinical training labels | Proteins | Spearman | Δ vs ours | AUC |
|---|---|---|---|---|---|
| **MissensExplorer v2.5** | ClinVar (other genes) | 61 | 0.462 [0.422–0.502] | — | 0.780 |
| AlphaMissense | none (population) | 59 | 0.453 [0.408–0.495] | -0.013 [-0.025, -0.001] | 0.776 |
| popEVE | none (population) | 58 | 0.453 [0.411–0.494] | -0.009 [-0.025, +0.005] | 0.767 |
| gMVP | ClinVar, HGMD, UniProt | 59 | 0.431 [0.387–0.469] | -0.035 [-0.053, -0.017] | 0.763 |
| MetaRNN | ClinVar | 60 | 0.431 [0.391–0.472] | -0.034 [-0.052, -0.018] | 0.762 |
| VARITY_R | ClinVar, HGMD | 57 | 0.430 [0.384–0.471] | -0.033 [-0.053, -0.015] | 0.759 |
| VARITY_R (leave-one-out) | ClinVar, HGMD | 57 | 0.425 [0.380–0.466] | -0.038 [-0.058, -0.019] | 0.756 |
| REVEL | HGMD | 61 | 0.419 [0.380–0.460] | -0.043 [-0.053, -0.032] | 0.755 |
| BayesDel (no AF) | ClinVar, UniProt | 61 | 0.411 [0.374–0.447] | -0.051 [-0.063, -0.039] | 0.746 |
| BayesDel (with AF) | ClinVar, UniProt | 61 | 0.410 [0.374–0.445] | -0.051 [-0.064, -0.039] | 0.746 |
| PHACTboost | ClinVar | 59 | 0.410 [0.372–0.449] | -0.050 [-0.069, -0.033] | 0.743 |
| EVE | none | 30 | 0.408 [0.350–0.467] | -0.039 [-0.072, -0.003] | 0.747 |
| MutFormer | HGMD | 61 | 0.408 [0.362–0.449] | -0.053 [-0.071, -0.037] | 0.745 |
| MutationAssessor | none | 57 | 0.408 [0.359–0.453] | -0.054 [-0.078, -0.033] | 0.744 |
| DEOGEN2 | UniProt | 59 | 0.405 [0.360–0.448] | -0.061 [-0.086, -0.037] | 0.744 |
| ClinPred | ClinVar | 61 | 0.402 [0.366–0.440] | -0.059 [-0.073, -0.044] | 0.744 |
| MutPred2 | HGMD, SwissVar | 61 | 0.402 [0.365–0.438] | -0.059 [-0.073, -0.045] | 0.747 |
| ESM1b | none | 61 | 0.400 [0.357–0.439] | -0.062 [-0.082, -0.042] | 0.740 |
| Eigen | none | 58 | 0.386 [0.345–0.423] | -0.076 [-0.096, -0.058] | 0.734 |
| MetaLR | UniProt | 61 | 0.369 [0.328–0.410] | -0.093 [-0.118, -0.066] | 0.723 |
| SIFT4G | none | 59 | 0.369 [0.337–0.403] | -0.095 [-0.116, -0.072] | 0.718 |
| PROVEAN | none | 60 | 0.359 [0.322–0.397] | -0.106 [-0.134, -0.080] | 0.714 |
| MPC | ClinVar (402 variants) | 58 | 0.352 [0.310–0.392] | -0.116 [-0.143, -0.091] | 0.715 |
| MVP | ClinVar, HGMD, UniProt | 61 | 0.340 [0.298–0.382] | -0.122 [-0.147, -0.096] | 0.703 |
| MetaSVM | UniProt | 61 | 0.339 [0.289–0.387] | -0.123 [-0.164, -0.084] | 0.710 |
| CADD | none (simulated) | 61 | 0.337 [0.298–0.373] | -0.125 [-0.152, -0.100] | 0.700 |
| MutationTaster2021 | ClinVar, HGMD (+ look-up) | 58 | 0.326 [0.285–0.366] | -0.136 [-0.161, -0.109] | 0.695 |
| PolyPhen-2 | UniProt | 50 | 0.325 [0.280–0.370] | -0.130 [-0.167, -0.097] | 0.685 |
| SIFT | none | 60 | 0.320 [0.280–0.360] | -0.141 [-0.172, -0.112] | 0.672 |
| FATHMM-XF | HGMD | 58 | 0.282 [0.246–0.321] | -0.180 [-0.208, -0.152] | 0.673 |
| PrimateAI | none (population) | 61 | 0.237 [0.201–0.276] | -0.225 [-0.254, -0.197] | 0.642 |
| LIST-S2 | few parameters | 60 | 0.236 [0.198–0.271] | -0.226 [-0.256, -0.195] | 0.636 |
| phyloP (100 vertebrates) | none | 61 | 0.234 [0.198–0.271] | -0.228 [-0.256, -0.199] | 0.642 |
| DANN | none (simulated) | 61 | 0.193 [0.157–0.230] | -0.269 [-0.302, -0.234] | 0.610 |

### Activity

| Tool | Proteins | Spearman | Δ vs ours |
|---|---|---|---|
| AlphaMissense | 14 | 0.531 | +0.009 [-0.010, +0.031] |
| MissensExplorer v2.5 | 14 | 0.523 | — |
| popEVE | 14 | 0.512 | -0.010 [-0.031, +0.006] |
| VARITY_R | 14 | 0.491 | -0.032 [-0.043, -0.021] |
| VARITY_R (leave-one-out) | 14 | 0.490 | -0.032 [-0.044, -0.021] |
| MutFormer | 14 | 0.481 | -0.042 [-0.056, -0.027] |
| ESM1b | 14 | 0.474 | -0.049 [-0.068, -0.029] |
| EVE | 11 | 0.474 | -0.052 [-0.096, -0.013] |
| MetaRNN | 14 | 0.470 | -0.053 [-0.082, -0.029] |
| PHACTboost | 14 | 0.465 | -0.058 [-0.091, -0.028] |
| gMVP | 14 | 0.464 | -0.058 [-0.085, -0.030] |
| REVEL | 14 | 0.464 | -0.059 [-0.076, -0.041] |
| MutationAssessor | 14 | 0.461 | -0.061 [-0.085, -0.034] |
| DEOGEN2 | 14 | 0.457 | -0.066 [-0.097, -0.036] |
| ClinPred | 14 | 0.455 | -0.067 [-0.095, -0.043] |
| MutPred2 | 14 | 0.450 | -0.073 [-0.100, -0.049] |
| BayesDel (with AF) | 14 | 0.446 | -0.077 [-0.096, -0.059] |
| BayesDel (no AF) | 14 | 0.443 | -0.080 [-0.098, -0.059] |
| Eigen | 13 | 0.436 | -0.078 [-0.102, -0.054] |
| CADD | 14 | 0.429 | -0.094 [-0.116, -0.073] |
| MPC | 13 | 0.417 | -0.115 [-0.152, -0.085] |
| PROVEAN | 14 | 0.413 | -0.110 [-0.139, -0.081] |
| SIFT4G | 14 | 0.406 | -0.117 [-0.150, -0.086] |
| PolyPhen-2 | 12 | 0.398 | -0.139 [-0.188, -0.088] |
| MVP | 14 | 0.389 | -0.134 [-0.175, -0.093] |
| MetaLR | 14 | 0.386 | -0.136 [-0.188, -0.087] |
| MutationTaster2021 | 14 | 0.357 | -0.158 [-0.185, -0.131] |
| SIFT | 14 | 0.356 | -0.166 [-0.233, -0.116] |
| MetaSVM | 14 | 0.348 | -0.174 [-0.246, -0.112] |
| FATHMM-XF | 13 | 0.328 | -0.186 [-0.227, -0.147] |
| PrimateAI | 14 | 0.314 | -0.209 [-0.264, -0.156] |
| LIST-S2 | 14 | 0.309 | -0.214 [-0.264, -0.172] |
| phyloP (100 vertebrates) | 14 | 0.283 | -0.240 [-0.289, -0.195] |
| DANN | 14 | 0.256 | -0.267 [-0.310, -0.225] |

### Binding

| Tool | Proteins | Spearman | Δ vs ours |
|---|---|---|---|
| DEOGEN2 | 6 | 0.458 | +0.002 [-0.040, +0.048] |
| MissensExplorer v2.5 | 6 | 0.455 | — |
| popEVE | 6 | 0.443 | -0.012 [-0.028, +0.008] |
| VARITY_R | 6 | 0.438 | -0.017 [-0.072, +0.027] |
| VARITY_R (leave-one-out) | 6 | 0.438 | -0.018 [-0.074, +0.028] |
| AlphaMissense | 6 | 0.436 | -0.020 [-0.062, +0.008] |
| MutationAssessor | 5 | 0.432 | -0.008 [-0.028, +0.010] |
| BayesDel (no AF) | 6 | 0.432 | -0.024 [-0.058, +0.015] |
| gMVP | 6 | 0.425 | -0.030 [-0.091, +0.022] |
| ClinPred | 6 | 0.419 | -0.036 [-0.067, -0.011] |
| REVEL | 6 | 0.419 | -0.036 [-0.077, -0.002] |
| BayesDel (with AF) | 6 | 0.415 | -0.040 [-0.083, +0.008] |
| MutPred2 | 6 | 0.413 | -0.043 [-0.116, +0.021] |
| PROVEAN | 6 | 0.407 | -0.048 [-0.081, -0.020] |
| MetaSVM | 6 | 0.399 | -0.056 [-0.134, +0.018] |
| MutFormer | 6 | 0.397 | -0.058 [-0.120, +0.002] |
| MetaRNN | 6 | 0.388 | -0.067 [-0.194, +0.019] |
| Eigen | 5 | 0.377 | -0.123 [-0.258, -0.038] |
| PHACTboost | 6 | 0.377 | -0.079 [-0.140, -0.016] |
| MPC | 6 | 0.363 | -0.093 [-0.171, -0.005] |
| SIFT4G | 6 | 0.346 | -0.109 [-0.194, -0.051] |
| MetaLR | 6 | 0.341 | -0.114 [-0.241, +0.016] |
| SIFT | 6 | 0.339 | -0.116 [-0.206, -0.031] |
| ESM1b | 6 | 0.329 | -0.126 [-0.242, -0.039] |
| CADD | 6 | 0.318 | -0.137 [-0.232, -0.048] |
| PolyPhen-2 | 5 | 0.316 | -0.134 [-0.349, -0.005] |
| MVP | 6 | 0.309 | -0.146 [-0.248, -0.057] |
| MutationTaster2021 | 6 | 0.293 | -0.167 [-0.251, -0.080] |
| FATHMM-XF | 5 | 0.293 | -0.208 [-0.319, -0.097] |
| phyloP (100 vertebrates) | 6 | 0.274 | -0.181 [-0.282, -0.080] |
| LIST-S2 | 6 | 0.266 | -0.190 [-0.286, -0.104] |
| PrimateAI | 6 | 0.186 | -0.269 [-0.374, -0.156] |
| DANN | 6 | 0.164 | -0.291 [-0.477, -0.124] |

### Expression

| Tool | Proteins | Spearman | Δ vs ours |
|---|---|---|---|
| popEVE | 12 | 0.475 | +0.008 [-0.015, +0.029] |
| MissensExplorer v2.5 | 12 | 0.467 | — |
| AlphaMissense | 12 | 0.459 | -0.012 [-0.036, +0.009] |
| MetaRNN | 11 | 0.457 | -0.028 [-0.054, -0.001] |
| VARITY_R | 11 | 0.442 | -0.048 [-0.116, -0.001] |
| ESM1b | 12 | 0.438 | -0.029 [-0.064, +0.010] |
| PHACTboost | 12 | 0.437 | -0.034 [-0.058, -0.010] |
| VARITY_R (leave-one-out) | 11 | 0.427 | -0.063 [-0.134, -0.009] |
| DEOGEN2 | 12 | 0.425 | -0.045 [-0.092, -0.002] |
| BayesDel (with AF) | 12 | 0.420 | -0.047 [-0.075, -0.007] |
| BayesDel (no AF) | 12 | 0.418 | -0.049 [-0.076, -0.012] |
| EVE | 6 | 0.416 | +0.008 [-0.049, +0.062] |
| REVEL | 12 | 0.413 | -0.054 [-0.071, -0.034] |
| MutPred2 | 12 | 0.408 | -0.059 [-0.078, -0.038] |
| PROVEAN | 11 | 0.408 | -0.078 [-0.115, -0.041] |
| MutFormer | 12 | 0.406 | -0.061 [-0.097, -0.027] |
| gMVP | 12 | 0.391 | -0.079 [-0.140, -0.032] |
| ClinPred | 12 | 0.383 | -0.084 [-0.108, -0.061] |
| SIFT4G | 11 | 0.377 | -0.108 [-0.147, -0.068] |
| MutationAssessor | 11 | 0.377 | -0.078 [-0.141, -0.023] |
| SIFT | 12 | 0.360 | -0.107 [-0.131, -0.085] |
| CADD | 12 | 0.354 | -0.113 [-0.164, -0.069] |
| MetaLR | 12 | 0.354 | -0.113 [-0.157, -0.061] |
| PolyPhen-2 | 11 | 0.352 | -0.152 [-0.176, -0.127] |
| Eigen | 12 | 0.350 | -0.117 [-0.158, -0.075] |
| MPC | 12 | 0.346 | -0.124 [-0.169, -0.079] |
| MutationTaster2021 | 11 | 0.329 | -0.152 [-0.219, -0.091] |
| FATHMM-XF | 12 | 0.321 | -0.146 [-0.199, -0.096] |
| MetaSVM | 12 | 0.321 | -0.146 [-0.230, -0.071] |
| MVP | 12 | 0.316 | -0.151 [-0.207, -0.086] |
| phyloP (100 vertebrates) | 12 | 0.240 | -0.227 [-0.282, -0.176] |
| LIST-S2 | 12 | 0.238 | -0.229 [-0.298, -0.161] |
| PrimateAI | 12 | 0.231 | -0.239 [-0.283, -0.195] |
| DANN | 12 | 0.202 | -0.265 [-0.348, -0.198] |

### OrganismalFitness

| Tool | Proteins | Spearman | Δ vs ours |
|---|---|---|---|
| MissensExplorer v2.5 | 22 | 0.416 | — |
| AlphaMissense | 21 | 0.413 | -0.014 [-0.035, +0.007] |
| VARITY_R | 21 | 0.395 | -0.017 [-0.045, +0.014] |
| gMVP | 21 | 0.394 | -0.032 [-0.062, -0.002] |
| popEVE | 21 | 0.394 | -0.018 [-0.056, +0.011] |
| VARITY_R (leave-one-out) | 21 | 0.388 | -0.024 [-0.054, +0.006] |
| MetaRNN | 22 | 0.386 | -0.030 [-0.047, -0.010] |
| ClinPred | 22 | 0.378 | -0.038 [-0.062, -0.007] |
| REVEL | 22 | 0.376 | -0.041 [-0.062, -0.019] |
| ESM1b | 22 | 0.373 | -0.043 [-0.061, -0.025] |
| MutPred2 | 22 | 0.371 | -0.045 [-0.064, -0.026] |
| EVE | 12 | 0.370 | -0.031 [-0.097, +0.040] |
| PHACTboost | 21 | 0.367 | -0.045 [-0.078, -0.014] |
| MutationAssessor | 21 | 0.366 | -0.061 [-0.093, -0.026] |
| Eigen | 21 | 0.362 | -0.051 [-0.075, -0.026] |
| BayesDel (with AF) | 22 | 0.361 | -0.055 [-0.077, -0.034] |
| BayesDel (no AF) | 22 | 0.360 | -0.056 [-0.078, -0.034] |
| SIFT4G | 22 | 0.357 | -0.059 [-0.097, -0.014] |
| MutFormer | 22 | 0.349 | -0.067 [-0.104, -0.032] |
| DEOGEN2 | 21 | 0.338 | -0.089 [-0.136, -0.042] |
| MutationTaster2021 | 21 | 0.330 | -0.070 [-0.112, -0.025] |
| MetaLR | 22 | 0.327 | -0.089 [-0.127, -0.054] |
| CADD | 22 | 0.322 | -0.094 [-0.140, -0.052] |
| MVP | 22 | 0.321 | -0.095 [-0.131, -0.058] |
| PROVEAN | 22 | 0.315 | -0.101 [-0.160, -0.042] |
| MPC | 22 | 0.299 | -0.117 [-0.181, -0.067] |
| PolyPhen-2 | 17 | 0.295 | -0.092 [-0.150, -0.047] |
| SIFT | 22 | 0.290 | -0.126 [-0.168, -0.085] |
| MetaSVM | 22 | 0.284 | -0.132 [-0.222, -0.063] |
| FATHMM-XF | 21 | 0.246 | -0.168 [-0.214, -0.121] |
| LIST-S2 | 21 | 0.232 | -0.182 [-0.238, -0.124] |
| PrimateAI | 22 | 0.226 | -0.190 [-0.246, -0.139] |
| phyloP (100 vertebrates) | 22 | 0.209 | -0.207 [-0.260, -0.153] |
| DANN | 22 | 0.184 | -0.232 [-0.285, -0.172] |

### Stability

| Tool | Proteins | Spearman | Δ vs ours |
|---|---|---|---|
| popEVE | 12 | 0.488 | +0.004 [-0.024, +0.030] |
| MissensExplorer v2.5 | 14 | 0.480 | — |
| gMVP | 13 | 0.477 | -0.000 [-0.021, +0.023] |
| MetaRNN | 14 | 0.464 | -0.017 [-0.053, +0.012] |
| MutationAssessor | 13 | 0.455 | -0.023 [-0.091, +0.041] |
| AlphaMissense | 13 | 0.449 | -0.028 [-0.049, -0.007] |
| MutFormer | 14 | 0.448 | -0.033 [-0.050, -0.017] |
| REVEL | 14 | 0.444 | -0.036 [-0.058, -0.016] |
| BayesDel (no AF) | 14 | 0.444 | -0.037 [-0.064, -0.008] |
| BayesDel (with AF) | 14 | 0.444 | -0.037 [-0.064, -0.006] |
| DEOGEN2 | 13 | 0.441 | -0.037 [-0.092, +0.016] |
| MetaLR | 14 | 0.429 | -0.051 [-0.103, -0.004] |
| Eigen | 14 | 0.420 | -0.061 [-0.084, -0.035] |
| VARITY_R (leave-one-out) | 12 | 0.414 | -0.059 [-0.095, -0.020] |
| VARITY_R | 12 | 0.413 | -0.060 [-0.100, -0.022] |
| PHACTboost | 13 | 0.413 | -0.064 [-0.110, -0.025] |
| MutPred2 | 14 | 0.405 | -0.075 [-0.113, -0.037] |
| MetaSVM | 14 | 0.402 | -0.078 [-0.131, -0.026] |
| ClinPred | 14 | 0.400 | -0.080 [-0.102, -0.059] |
| ESM1b | 14 | 0.383 | -0.098 [-0.146, -0.048] |
| MPC | 12 | 0.364 | -0.143 [-0.189, -0.105] |
| SIFT4G | 13 | 0.361 | -0.116 [-0.152, -0.078] |
| MVP | 14 | 0.352 | -0.128 [-0.192, -0.070] |
| PROVEAN | 14 | 0.326 | -0.156 [-0.203, -0.113] |
| SIFT | 13 | 0.297 | -0.180 [-0.263, -0.095] |
| CADD | 14 | 0.285 | -0.196 [-0.250, -0.137] |
| MutationTaster2021 | 13 | 0.279 | -0.217 [-0.257, -0.178] |
| PolyPhen-2 | 11 | 0.276 | -0.180 [-0.279, -0.093] |
| FATHMM-XF | 14 | 0.271 | -0.210 [-0.272, -0.149] |
| PrimateAI | 14 | 0.222 | -0.259 [-0.311, -0.199] |
| phyloP (100 vertebrates) | 14 | 0.212 | -0.268 [-0.326, -0.211] |
| LIST-S2 | 14 | 0.182 | -0.299 [-0.352, -0.247] |
| DANN | 14 | 0.176 | -0.305 [-0.353, -0.250] |

Assays left out (25): their sequence does not occur exactly in the gene's canonical protein (isolated domains with tags, engineered constructs) — AMFR_HUMAN_Tsuboyama_2023_4G3O, CASP3_HUMAN_Roychowdhury_2020, CASP7_HUMAN_Roychowdhury_2020, CCR5_HUMAN_Gill_2023, GLPA_HUMAN_Elazar_2016, HECD1_HUMAN_Tsuboyama_2023_3DKM, ILF3_HUMAN_Tsuboyama_2023_2L33, KCNE1_HUMAN_Muhammad_2023_expression, KCNE1_HUMAN_Muhammad_2023_function, NKX31_HUMAN_Tsuboyama_2023_2L9R, P53_HUMAN_Giacomelli_2018_Null_Etoposide, P53_HUMAN_Giacomelli_2018_Null_Nutlin, P53_HUMAN_Giacomelli_2018_WT_Nutlin, PITX2_HUMAN_Tsuboyama_2023_2L7M, PKN1_HUMAN_Tsuboyama_2023_1URF, PPARG_HUMAN_Majithia_2016, Q53Z42_HUMAN_McShan_2019_binding-TAPBPR, Q53Z42_HUMAN_McShan_2019_expression, RBP1_HUMAN_Tsuboyama_2023_2KWH, SCN5A_HUMAN_Glazer_2019, SRBS1_HUMAN_Tsuboyama_2023_2O2W, SYUA_HUMAN_Newberry_2020, TNKS2_HUMAN_Tsuboyama_2023_5JRT, TPOR_HUMAN_Bridgford_2020, UBC9_HUMAN_Weile_2017.
