# Validation on deep mutational scanning data (ProteinGym)

13408 single-nucleotide missense variants from 71 human assays on 61 proteins of [ProteinGym](https://marks.hms.harvard.edu/proteingym/ProteinGym_v1.3/DMS_ProteinGym_substitutions.zip) (up to 200 per assay). No ClinVar label is involved: each variant's effect was measured in the laboratory. Produced by `benchmarks/validate_proteingym.py`.

- **Spearman**: correlation between the tool's score and the measured loss of function, averaged over the assays of each protein and then over proteins. 95% bootstrap interval over proteins.
- **Δ vs ours**: the tool's mean minus MissensExplorer's on the same proteins and variants; an interval that excludes 0 is a real difference.
- **AUC**: against ProteinGym's binary label (fit/unfit), same averaging.
- DMS assays measure one molecular function in one experimental system, not disease, so every tool's correlation is moderate. VARITY used some DMS data in training.

| Tool | Clinical training labels | Proteins | Spearman | Δ vs ours | AUC |
|---|---|---|---|---|---|
| **MissensExplorer v2.3** | ClinVar (other genes) | 61 | 0.469 [0.430–0.507] | — | 0.783 |
| popEVE | none (population) | 58 | 0.460 [0.420–0.500] | -0.009 [-0.029, +0.008] | 0.770 |
| AlphaMissense | none (population) | 59 | 0.458 [0.412–0.499] | -0.012 [-0.025, +0.000] | 0.776 |
| MetaRNN | ClinVar | 60 | 0.438 [0.398–0.477] | -0.035 [-0.052, -0.018] | 0.765 |
| VARITY_R | ClinVar, HGMD | 57 | 0.431 [0.386–0.472] | -0.037 [-0.059, -0.018] | 0.760 |
| gMVP | ClinVar, HGMD, UniProt | 59 | 0.430 [0.386–0.468] | -0.039 [-0.059, -0.021] | 0.761 |
| REVEL | HGMD | 61 | 0.425 [0.387–0.463] | -0.044 [-0.056, -0.034] | 0.759 |
| VARITY_R (leave-one-out) | ClinVar, HGMD | 57 | 0.424 [0.380–0.465] | -0.043 [-0.063, -0.024] | 0.756 |
| BayesDel (no AF) | ClinVar, UniProt | 61 | 0.412 [0.374–0.450] | -0.057 [-0.070, -0.044] | 0.748 |
| EVE | none | 30 | 0.412 [0.352–0.469] | -0.033 [-0.074, +0.006] | 0.747 |
| BayesDel (with AF) | ClinVar, UniProt | 61 | 0.411 [0.374–0.448] | -0.058 [-0.071, -0.045] | 0.748 |
| PHACTboost | ClinVar | 59 | 0.411 [0.375–0.448] | -0.054 [-0.076, -0.035] | 0.744 |
| MutationAssessor | none | 57 | 0.410 [0.361–0.457] | -0.058 [-0.084, -0.033] | 0.744 |
| DEOGEN2 | UniProt | 59 | 0.410 [0.366–0.452] | -0.062 [-0.088, -0.036] | 0.745 |
| MutFormer | HGMD | 61 | 0.410 [0.365–0.451] | -0.059 [-0.078, -0.041] | 0.748 |
| ESM1b | none | 61 | 0.408 [0.368–0.449] | -0.061 [-0.083, -0.039] | 0.745 |
| MutPred2 | HGMD, SwissVar | 61 | 0.406 [0.368–0.441] | -0.063 [-0.078, -0.049] | 0.747 |
| ClinPred | ClinVar | 61 | 0.405 [0.367–0.442] | -0.065 [-0.080, -0.048] | 0.745 |
| Eigen | none | 58 | 0.394 [0.356–0.432] | -0.077 [-0.097, -0.058] | 0.737 |
| MetaLR | UniProt | 61 | 0.373 [0.330–0.414] | -0.096 [-0.122, -0.069] | 0.726 |
| SIFT4G | none | 59 | 0.372 [0.342–0.404] | -0.097 [-0.120, -0.072] | 0.719 |
| PROVEAN | none | 59 | 0.367 [0.330–0.402] | -0.103 [-0.130, -0.074] | 0.715 |
| MPC | ClinVar (402 variants) | 58 | 0.359 [0.317–0.400] | -0.117 [-0.145, -0.090] | 0.719 |
| CADD | none (simulated) | 61 | 0.344 [0.304–0.382] | -0.126 [-0.154, -0.100] | 0.702 |
| MetaSVM | UniProt | 61 | 0.342 [0.291–0.391] | -0.127 [-0.168, -0.088] | 0.713 |
| MVP | ClinVar, HGMD, UniProt | 61 | 0.342 [0.300–0.383] | -0.127 [-0.155, -0.099] | 0.708 |
| MutationTaster2021 | ClinVar, HGMD (+ look-up) | 55 | 0.334 [0.297–0.370] | -0.133 [-0.160, -0.106] | 0.700 |
| PolyPhen-2 | UniProt | 50 | 0.330 [0.286–0.377] | -0.135 [-0.176, -0.098] | 0.686 |
| SIFT | none | 60 | 0.323 [0.284–0.363] | -0.145 [-0.177, -0.116] | 0.671 |
| FATHMM-XF | HGMD | 58 | 0.286 [0.246–0.328] | -0.185 [-0.212, -0.156] | 0.674 |
| LIST-S2 | few parameters | 60 | 0.246 [0.204–0.285] | -0.223 [-0.256, -0.189] | 0.642 |
| phyloP (100 vertebrates) | none | 61 | 0.237 [0.197–0.276] | -0.232 [-0.262, -0.202] | 0.644 |
| PrimateAI | none (population) | 61 | 0.233 [0.195–0.268] | -0.236 [-0.266, -0.206] | 0.638 |
| DANN | none (simulated) | 61 | 0.197 [0.159–0.234] | -0.272 [-0.310, -0.233] | 0.610 |

Assays left out (25): their sequence does not occur exactly in the gene's canonical protein (isolated domains with tags, engineered constructs) — AMFR_HUMAN_Tsuboyama_2023_4G3O, CASP3_HUMAN_Roychowdhury_2020, CASP7_HUMAN_Roychowdhury_2020, CCR5_HUMAN_Gill_2023, GLPA_HUMAN_Elazar_2016, HECD1_HUMAN_Tsuboyama_2023_3DKM, ILF3_HUMAN_Tsuboyama_2023_2L33, KCNE1_HUMAN_Muhammad_2023_expression, KCNE1_HUMAN_Muhammad_2023_function, NKX31_HUMAN_Tsuboyama_2023_2L9R, P53_HUMAN_Giacomelli_2018_Null_Etoposide, P53_HUMAN_Giacomelli_2018_Null_Nutlin, P53_HUMAN_Giacomelli_2018_WT_Nutlin, PITX2_HUMAN_Tsuboyama_2023_2L7M, PKN1_HUMAN_Tsuboyama_2023_1URF, PPARG_HUMAN_Majithia_2016, Q53Z42_HUMAN_McShan_2019_binding-TAPBPR, Q53Z42_HUMAN_McShan_2019_expression, RBP1_HUMAN_Tsuboyama_2023_2KWH, SCN5A_HUMAN_Glazer_2019, SRBS1_HUMAN_Tsuboyama_2023_2O2W, SYUA_HUMAN_Newberry_2020, TNKS2_HUMAN_Tsuboyama_2023_5JRT, TPOR_HUMAN_Bridgford_2020, UBC9_HUMAN_Weile_2017.
