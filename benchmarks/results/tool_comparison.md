# MissensExplorer and other missense predictors

Same held-out benchmark as [clinvar_benchmark.md](clinvar_benchmark.md): 1972 expert-reviewed ClinVar missense variants from 1181 genes, none of them used to train the MissensExplorer model. Other tools' scores come from dbNSFP through the Ensembl VEP REST API (Ensembl VEP REST, release 116). Produced by `benchmarks/compare_tools.py`.

- **AUC** (area under the ROC curve) is threshold-free: the chance that a random pathogenic variant scores higher than a random benign one. 95% bootstrap intervals in brackets.
- **Δ vs ours**: the tool's AUC minus MissensExplorer's on the same variants, with a paired bootstrap interval. An interval that excludes 0 is a real difference on this benchmark.
- **Rare**: gnomAD grpmax < 0.1% (964 pathogenic, 199 benign). **New since 2023**: variants whose pathogenic/benign label was not yet in the January 2023 ClinVar release (282 pathogenic, 172 benign); tools trained before 2023 cannot have seen these labels.
- **Clinical training labels**: the disease-variant databases a tool was trained on. A tool trained on ClinVar may have seen some of these very variants (all such tools were trained on ClinVar releases from 2020 or earlier); HGMD and UniProt overlap ClinVar's pathogenic variants. *none* = no clinical labels.
- **Uses allele frequency**: the tool takes gnomAD/1000 Genomes frequencies as input. ClinVar curators use the same frequencies to call variants benign (ACMG BA1/BS1), so these tools look better on ClinVar than the protein-level evidence alone would. MissensExplorer applies frequency as a separate, visible rule instead; the line *+ frequency rule (BA1)* scores the default verdict, and *gnomAD frequency alone* shows how far frequency goes by itself.
- **MutationTaster2021** calls variants found in ClinVar or HGMD disease-causing by look-up, so on this benchmark it partly reads the answer.

## All variants

| Tool | Clinical training labels | Uses allele frequency | Coverage | AUC | Δ vs ours |
|---|---|---|---|---|---|
| MetaRNN | ClinVar | yes | 99% | 0.996 [0.994–0.998] | +0.011 [+0.007, +0.015] |
| **MissensExplorer v2.3 + frequency as a feature (exploratory)** | ClinVar (other genes) | yes | 100% | 0.996 [0.994–0.998] | +0.010 [+0.007, +0.014] |
| ClinPred | ClinVar | yes | 99% | 0.993 [0.989–0.996] | +0.008 [+0.004, +0.011] |
| BayesDel (with AF) | ClinVar, UniProt | yes | 100% | 0.992 [0.988–0.995] | +0.007 [+0.003, +0.011] |
| MutationTaster2021 | ClinVar, HGMD (+ look-up) | no | 96% | 0.992 [0.988–0.996] | +0.006 [+0.002, +0.011] |
| **MissensExplorer v2.3 + frequency rule (BA1)** | ClinVar (other genes) | yes | 100% | 0.986 [0.982–0.990] | +0.001 [+0.000, +0.001] |
| **MissensExplorer v2.3** | ClinVar (other genes) | no | 100% | 0.985 [0.980–0.990] | — |
| PHACTboost | ClinVar | no | 93% | 0.983 [0.977–0.988] | -0.003 [-0.007, +0.001] |
| REVEL | HGMD | no | 96% | 0.982 [0.977–0.986] | -0.004 [-0.007, -0.001] |
| gnomAD frequency alone (reference) | — | yes | 100% | 0.978 [0.973–0.983] | -0.007 [-0.013, -0.001] |
| BayesDel (no AF) | ClinVar, UniProt | no | 100% | 0.977 [0.971–0.983] | -0.008 [-0.013, -0.004] |
| gMVP | ClinVar, HGMD, UniProt | no | 93% | 0.972 [0.964–0.978] | -0.015 [-0.021, -0.009] |
| VARITY_R | ClinVar, HGMD | no | 91% | 0.971 [0.964–0.977] | -0.014 [-0.019, -0.009] |
| VARITY_R (leave-one-out) | ClinVar, HGMD | no | 91% | 0.968 [0.961–0.976] | -0.017 [-0.022, -0.012] |
| AlphaMissense | none (population) | no | 99% | 0.964 [0.956–0.970] | -0.022 [-0.027, -0.018] |
| MutPred2 | HGMD, SwissVar | no | 100% | 0.958 [0.949–0.966] | -0.027 [-0.033, -0.021] |
| DEOGEN2 | UniProt | no | 97% | 0.954 [0.946–0.963] | -0.031 [-0.038, -0.023] |
| MetaSVM | UniProt | yes | 99% | 0.954 [0.945–0.963] | -0.031 [-0.039, -0.023] |
| MetaLR | UniProt | yes | 99% | 0.952 [0.943–0.961] | -0.033 [-0.041, -0.025] |
| MutFormer | HGMD | no | 99% | 0.952 [0.942–0.960] | -0.034 [-0.043, -0.025] |
| MVP | ClinVar, HGMD, UniProt | no | 89% | 0.952 [0.941–0.962] | -0.032 [-0.041, -0.022] |
| ESM1b | none | no | 100% | 0.935 [0.925–0.946] | -0.050 [-0.059, -0.040] |
| Eigen | none | no | 92% | 0.933 [0.922–0.943] | -0.053 [-0.063, -0.045] |
| EVE | none | no | 59% | 0.928 [0.913–0.942] | -0.054 [-0.068, -0.042] |
| CADD | none (simulated) | no | 100% | 0.923 [0.912–0.935] | -0.062 [-0.072, -0.052] |
| FATHMM-XF | HGMD | no | 92% | 0.914 [0.901–0.926] | -0.072 [-0.083, -0.061] |
| MutationAssessor | none | no | 92% | 0.912 [0.899–0.925] | -0.073 [-0.085, -0.061] |
| PROVEAN | none | no | 95% | 0.910 [0.896–0.923] | -0.076 [-0.088, -0.064] |
| popEVE | none (population) | no | 83% | 0.906 [0.892–0.920] | -0.079 [-0.091, -0.066] |
| SIFT4G | none | no | 95% | 0.898 [0.884–0.912] | -0.087 [-0.100, -0.075] |
| LIST-S2 | few parameters | no | 98% | 0.893 [0.879–0.907] | -0.092 [-0.105, -0.079] |
| PolyPhen-2 | UniProt | no | 93% | 0.884 [0.867–0.898] | -0.101 [-0.118, -0.087] |
| PrimateAI | none (population) | no | 98% | 0.877 [0.863–0.892] | -0.108 [-0.122, -0.094] |
| SIFT | none | no | 98% | 0.862 [0.845–0.877] | -0.123 [-0.138, -0.110] |
| phyloP (100 vertebrates) | none | no | 100% | 0.859 [0.843–0.874] | -0.126 [-0.142, -0.112] |
| DANN | none (simulated) | no | 100% | 0.827 [0.808–0.846] | -0.158 [-0.176, -0.140] |
| MPC | ClinVar (402 variants) | no | 86% | 0.817 [0.797–0.837] | -0.169 [-0.189, -0.150] |

## Rare variants

| Tool | Clinical training labels | Uses allele frequency | Coverage | AUC | Δ vs ours |
|---|---|---|---|---|---|
| MetaRNN | ClinVar | yes | 99% | 0.983 [0.973–0.991] | +0.009 [+0.002, +0.018] |
| **MissensExplorer v2.3 + frequency as a feature (exploratory)** | ClinVar (other genes) | yes | 100% | 0.983 [0.973–0.990] | +0.009 [+0.003, +0.015] |
| **MissensExplorer v2.3** | ClinVar (other genes) | no | 100% | 0.974 [0.961–0.985] | — |
| **MissensExplorer v2.3 + frequency rule (BA1)** | ClinVar (other genes) | yes | 100% | 0.974 [0.961–0.985] | +0.000 [+0.000, +0.000] |
| MutationTaster2021 | ClinVar, HGMD (+ look-up) | no | 96% | 0.974 [0.956–0.988] | -0.002 [-0.016, +0.012] |
| ClinPred | ClinVar | yes | 99% | 0.973 [0.958–0.985] | -0.002 [-0.010, +0.006] |
| BayesDel (with AF) | ClinVar, UniProt | yes | 100% | 0.972 [0.955–0.986] | -0.002 [-0.014, +0.007] |
| REVEL | HGMD | no | 96% | 0.972 [0.959–0.983] | -0.004 [-0.009, +0.001] |
| gMVP | ClinVar, HGMD, UniProt | no | 93% | 0.966 [0.951–0.979] | -0.008 [-0.019, +0.002] |
| BayesDel (no AF) | ClinVar, UniProt | no | 100% | 0.962 [0.943–0.977] | -0.013 [-0.025, -0.003] |
| VARITY_R | ClinVar, HGMD | no | 91% | 0.960 [0.944–0.974] | -0.013 [-0.022, -0.005] |
| VARITY_R (leave-one-out) | ClinVar, HGMD | no | 91% | 0.957 [0.940–0.971] | -0.017 [-0.025, -0.007] |
| PHACTboost | ClinVar | no | 93% | 0.957 [0.937–0.974] | -0.018 [-0.031, -0.004] |
| AlphaMissense | none (population) | no | 99% | 0.953 [0.938–0.966] | -0.022 [-0.032, -0.013] |
| MutPred2 | HGMD, SwissVar | no | 100% | 0.946 [0.927–0.962] | -0.028 [-0.043, -0.016] |
| MVP | ClinVar, HGMD, UniProt | no | 89% | 0.938 [0.917–0.958] | -0.037 [-0.055, -0.020] |
| DEOGEN2 | UniProt | no | 97% | 0.930 [0.910–0.950] | -0.044 [-0.060, -0.029] |
| MutationAssessor | none | no | 92% | 0.923 [0.904–0.943] | -0.050 [-0.065, -0.034] |
| ESM1b | none | no | 100% | 0.920 [0.899–0.942] | -0.054 [-0.072, -0.036] |
| Eigen | none | no | 92% | 0.916 [0.894–0.938] | -0.058 [-0.074, -0.041] |
| MetaSVM | UniProt | yes | 99% | 0.916 [0.887–0.940] | -0.058 [-0.083, -0.038] |
| gnomAD frequency alone (reference) | — | yes | 100% | 0.909 [0.890–0.928] | -0.065 [-0.085, -0.044] |
| MetaLR | UniProt | yes | 99% | 0.909 [0.882–0.931] | -0.065 [-0.087, -0.047] |
| PROVEAN | none | no | 95% | 0.905 [0.876–0.928] | -0.070 [-0.094, -0.049] |
| EVE | none | no | 59% | 0.901 [0.868–0.929] | -0.065 [-0.093, -0.040] |
| SIFT4G | none | no | 95% | 0.897 [0.872–0.920] | -0.077 [-0.103, -0.056] |
| FATHMM-XF | HGMD | no | 92% | 0.897 [0.869–0.921] | -0.077 [-0.101, -0.056] |
| MutFormer | HGMD | no | 99% | 0.889 [0.862–0.916] | -0.085 [-0.109, -0.060] |
| CADD | none (simulated) | no | 100% | 0.885 [0.855–0.912] | -0.089 [-0.116, -0.066] |
| PolyPhen-2 | UniProt | no | 93% | 0.883 [0.859–0.907] | -0.092 [-0.113, -0.070] |
| popEVE | none (population) | no | 83% | 0.856 [0.821–0.887] | -0.116 [-0.147, -0.089] |
| LIST-S2 | few parameters | no | 98% | 0.851 [0.822–0.879] | -0.123 [-0.151, -0.097] |
| phyloP (100 vertebrates) | none | no | 100% | 0.835 [0.802–0.864] | -0.139 [-0.168, -0.112] |
| SIFT | none | no | 98% | 0.805 [0.765–0.840] | -0.169 [-0.206, -0.136] |
| PrimateAI | none (population) | no | 98% | 0.793 [0.756–0.831] | -0.180 [-0.216, -0.147] |
| DANN | none (simulated) | no | 100% | 0.787 [0.748–0.826] | -0.187 [-0.224, -0.151] |
| MPC | ClinVar (402 variants) | no | 86% | 0.746 [0.705–0.783] | -0.229 [-0.268, -0.192] |

## New since 2023 variants

| Tool | Clinical training labels | Uses allele frequency | Coverage | AUC | Δ vs ours |
|---|---|---|---|---|---|
| MetaRNN | ClinVar | yes | 99% | 0.983 [0.971–0.992] | +0.012 [+0.002, +0.024] |
| **MissensExplorer v2.3 + frequency as a feature (exploratory)** | ClinVar (other genes) | yes | 100% | 0.982 [0.971–0.992] | +0.011 [+0.002, +0.022] |
| MutationTaster2021 | ClinVar, HGMD (+ look-up) | no | 96% | 0.972 [0.952–0.987] | +0.000 [-0.017, +0.016] |
| ClinPred | ClinVar | yes | 99% | 0.971 [0.953–0.986] | +0.001 [-0.011, +0.012] |
| **MissensExplorer v2.3 + frequency rule (BA1)** | ClinVar (other genes) | yes | 100% | 0.971 [0.952–0.986] | +0.000 [+0.000, +0.001] |
| **MissensExplorer v2.3** | ClinVar (other genes) | no | 100% | 0.971 [0.952–0.985] | — |
| BayesDel (with AF) | ClinVar, UniProt | yes | 100% | 0.969 [0.949–0.986] | -0.002 [-0.016, +0.010] |
| REVEL | HGMD | no | 96% | 0.966 [0.947–0.981] | -0.005 [-0.014, +0.003] |
| AlphaMissense | none (population) | no | 99% | 0.957 [0.940–0.974] | -0.014 [-0.026, -0.004] |
| PHACTboost | ClinVar | no | 93% | 0.956 [0.935–0.974] | -0.015 [-0.029, -0.004] |
| gMVP | ClinVar, HGMD, UniProt | no | 93% | 0.956 [0.934–0.975] | -0.018 [-0.032, -0.004] |
| VARITY_R | ClinVar, HGMD | no | 91% | 0.954 [0.932–0.972] | -0.016 [-0.031, -0.003] |
| VARITY_R (leave-one-out) | ClinVar, HGMD | no | 91% | 0.953 [0.931–0.973] | -0.017 [-0.031, -0.003] |
| BayesDel (no AF) | ClinVar, UniProt | no | 100% | 0.950 [0.926–0.970] | -0.021 [-0.037, -0.008] |
| MutPred2 | HGMD, SwissVar | no | 100% | 0.942 [0.919–0.964] | -0.029 [-0.046, -0.014] |
| gnomAD frequency alone (reference) | — | yes | 100% | 0.938 [0.919–0.957] | -0.033 [-0.057, -0.009] |
| MVP | ClinVar, HGMD, UniProt | no | 89% | 0.929 [0.899–0.954] | -0.042 [-0.068, -0.019] |
| MutationAssessor | none | no | 92% | 0.923 [0.897–0.948] | -0.046 [-0.066, -0.027] |
| Eigen | none | no | 92% | 0.921 [0.892–0.944] | -0.052 [-0.071, -0.033] |
| MutFormer | HGMD | no | 99% | 0.917 [0.889–0.942] | -0.054 [-0.079, -0.031] |
| DEOGEN2 | UniProt | no | 97% | 0.912 [0.884–0.940] | -0.058 [-0.082, -0.035] |
| FATHMM-XF | HGMD | no | 92% | 0.910 [0.879–0.937] | -0.063 [-0.088, -0.038] |
| ESM1b | none | no | 100% | 0.909 [0.881–0.938] | -0.062 [-0.088, -0.037] |
| MetaSVM | UniProt | yes | 99% | 0.907 [0.874–0.938] | -0.064 [-0.092, -0.037] |
| MetaLR | UniProt | yes | 99% | 0.906 [0.875–0.934] | -0.065 [-0.091, -0.041] |
| EVE | none | no | 59% | 0.906 [0.867–0.941] | -0.059 [-0.092, -0.029] |
| PolyPhen-2 | UniProt | no | 93% | 0.896 [0.864–0.928] | -0.077 [-0.107, -0.047] |
| PROVEAN | none | no | 95% | 0.893 [0.857–0.926] | -0.077 [-0.108, -0.049] |
| CADD | none (simulated) | no | 100% | 0.890 [0.859–0.920] | -0.081 [-0.108, -0.056] |
| SIFT4G | none | no | 95% | 0.886 [0.854–0.916] | -0.083 [-0.116, -0.054] |
| popEVE | none (population) | no | 83% | 0.867 [0.828–0.903] | -0.098 [-0.132, -0.068] |
| LIST-S2 | few parameters | no | 98% | 0.864 [0.828–0.896] | -0.107 [-0.139, -0.076] |
| phyloP (100 vertebrates) | none | no | 100% | 0.856 [0.821–0.890] | -0.115 [-0.145, -0.084] |
| SIFT | none | no | 98% | 0.812 [0.770–0.854] | -0.158 [-0.197, -0.119] |
| DANN | none (simulated) | no | 100% | 0.810 [0.767–0.850] | -0.161 [-0.201, -0.122] |
| PrimateAI | none (population) | no | 98% | 0.802 [0.753–0.842] | -0.169 [-0.214, -0.129] |
| MPC | ClinVar (402 variants) | no | 86% | 0.734 [0.680–0.784] | -0.236 [-0.286, -0.189] |

## Same variants for every tool (1290 variants scored by all tools with ≥ 90% coverage)

| Tool | All | Rare | New since 2023 |
|---|---|---|---|
| MetaRNN | 0.997 | 0.986 | 0.987 |
| MissensExplorer v2.3 + frequency as a feature (exploratory) | 0.997 | 0.982 | 0.984 |
| ClinPred | 0.994 | 0.973 | 0.978 |
| BayesDel (with AF) | 0.994 | 0.975 | 0.974 |
| MutationTaster2021 | 0.994 | 0.980 | 0.986 |
| MissensExplorer v2.3 + frequency rule (BA1) | 0.987 | 0.975 | 0.974 |
| MissensExplorer v2.3 | 0.987 | 0.975 | 0.974 |
| PHACTboost | 0.984 | 0.954 | 0.958 |
| REVEL | 0.984 | 0.973 | 0.971 |
| gnomAD frequency alone (reference) | 0.981 | 0.904 | 0.940 |
| BayesDel (no AF) | 0.980 | 0.964 | 0.954 |
| VARITY_R | 0.973 | 0.966 | 0.958 |
| VARITY_R (leave-one-out) | 0.970 | 0.961 | 0.957 |
| gMVP | 0.969 | 0.964 | 0.952 |
| AlphaMissense | 0.964 | 0.953 | 0.956 |
| MetaSVM | 0.960 | 0.919 | 0.926 |
| MutPred2 | 0.959 | 0.946 | 0.947 |
| MutFormer | 0.958 | 0.899 | 0.925 |
| DEOGEN2 | 0.956 | 0.934 | 0.926 |
| MetaLR | 0.955 | 0.909 | 0.923 |
| Eigen | 0.938 | 0.920 | 0.925 |
| ESM1b | 0.934 | 0.914 | 0.911 |
| CADD | 0.922 | 0.874 | 0.885 |
| PROVEAN | 0.914 | 0.915 | 0.911 |
| FATHMM-XF | 0.912 | 0.891 | 0.904 |
| MutationAssessor | 0.909 | 0.917 | 0.925 |
| SIFT4G | 0.907 | 0.899 | 0.904 |
| PolyPhen-2 | 0.896 | 0.889 | 0.907 |
| LIST-S2 | 0.891 | 0.845 | 0.868 |
| SIFT | 0.877 | 0.835 | 0.828 |
| PrimateAI | 0.876 | 0.804 | 0.808 |
| phyloP (100 vertebrates) | 0.863 | 0.842 | 0.863 |
| DANN | 0.822 | 0.749 | 0.790 |

## Calls at ClinGen-calibrated thresholds (PP3/BP4, supporting level), rare variants

| Tool | Pathogenic called damaging | Benign called benign | Pathogenic called benign | Benign called damaging | Decisive calls | Decisive calls right |
|---|---|---|---|---|---|---|
| **MissensExplorer v2.3 (default verdict)** | 91.6% | 79.4% | 2.3% | 6.0% | 92.4% | 96.8% |
| AlphaMissense | 71.5% | 68.8% | 3.4% | 4.0% | 74.5% | 95.3% |
| REVEL | 86.0% | 69.8% | 1.2% | 6.0% | 85.3% | 97.6% |
| CADD | 78.5% | 60.8% | 6.1% | 18.6% | 83.7% | 90.1% |
| SIFT | 80.5% | 41.7% | 2.4% | 27.6% | 80.6% | 91.7% |
| PolyPhen-2 | 56.5% | 49.2% | 6.3% | 5.5% | 61.5% | 89.9% |
| BayesDel (no AF) | 88.2% | 69.8% | 2.1% | 7.0% | 88.0% | 96.7% |
| VARITY_R (leave-one-out) | 77.7% | 60.8% | 2.1% | 5.0% | 77.4% | 96.7% |
| MutPred2 | 74.4% | 78.4% | 3.5% | 6.0% | 79.0% | 95.0% |
| MPC | 35.1% | 0.0% | 0.0% | 11.1% | 31.0% | 93.9% |
| PrimateAI | 43.5% | 46.2% | 6.0% | 14.1% | 51.3% | 85.6% |
| ESM1b | 74.3% | 67.8% | 5.3% | 7.5% | 78.8% | 92.8% |
| phyloP (100 vertebrates) | 54.5% | 42.7% | 6.5% | 9.5% | 59.5% | 88.2% |
