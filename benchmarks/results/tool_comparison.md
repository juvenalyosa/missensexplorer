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
| **MissensExplorer v2.5 + frequency as a feature (exploratory)** | ClinVar (other genes) | yes | 100% | 0.996 [0.994–0.997] | +0.010 [+0.007, +0.014] |
| ClinPred | ClinVar | yes | 99% | 0.993 [0.989–0.996] | +0.007 [+0.004, +0.011] |
| BayesDel (with AF) | ClinVar, UniProt | yes | 100% | 0.992 [0.988–0.995] | +0.007 [+0.003, +0.011] |
| MutationTaster2021 | ClinVar, HGMD (+ look-up) | no | 96% | 0.992 [0.988–0.995] | +0.006 [+0.002, +0.011] |
| **MissensExplorer v2.5 + frequency rule (BA1)** | ClinVar (other genes) | yes | 100% | 0.986 [0.982–0.990] | +0.001 [+0.000, +0.001] |
| **MissensExplorer v2.5** | ClinVar (other genes) | no | 100% | 0.985 [0.981–0.990] | — |
| PHACTboost | ClinVar | no | 93% | 0.983 [0.978–0.988] | -0.003 [-0.007, +0.001] |
| REVEL | HGMD | no | 96% | 0.982 [0.977–0.987] | -0.004 [-0.007, -0.001] |
| gnomAD frequency alone (reference) | — | yes | 100% | 0.978 [0.973–0.983] | -0.007 [-0.013, -0.001] |
| BayesDel (no AF) | ClinVar, UniProt | no | 100% | 0.977 [0.971–0.983] | -0.008 [-0.013, -0.004] |
| gMVP | ClinVar, HGMD, UniProt | no | 93% | 0.972 [0.965–0.978] | -0.015 [-0.021, -0.010] |
| VARITY_R | ClinVar, HGMD | no | 91% | 0.971 [0.964–0.977] | -0.014 [-0.019, -0.010] |
| VARITY_R (leave-one-out) | ClinVar, HGMD | no | 91% | 0.968 [0.961–0.975] | -0.017 [-0.022, -0.012] |
| AlphaMissense | none (population) | no | 99% | 0.964 [0.956–0.970] | -0.022 [-0.028, -0.017] |
| MutPred2 | HGMD, SwissVar | no | 100% | 0.958 [0.950–0.966] | -0.027 [-0.034, -0.021] |
| DEOGEN2 | UniProt | no | 97% | 0.954 [0.945–0.963] | -0.031 [-0.039, -0.023] |
| MetaSVM | UniProt | yes | 99% | 0.954 [0.944–0.963] | -0.032 [-0.040, -0.024] |
| MetaLR | UniProt | yes | 99% | 0.952 [0.944–0.961] | -0.033 [-0.041, -0.026] |
| MutFormer | HGMD | no | 99% | 0.952 [0.943–0.960] | -0.034 [-0.042, -0.025] |
| MVP | ClinVar, HGMD, UniProt | no | 89% | 0.952 [0.942–0.961] | -0.032 [-0.042, -0.023] |
| ESM1b | none | no | 100% | 0.935 [0.923–0.946] | -0.050 [-0.061, -0.040] |
| Eigen | none | no | 92% | 0.933 [0.922–0.944] | -0.054 [-0.063, -0.045] |
| EVE | none | no | 59% | 0.928 [0.911–0.942] | -0.054 [-0.069, -0.041] |
| CADD | none (simulated) | no | 100% | 0.923 [0.911–0.934] | -0.062 [-0.073, -0.052] |
| FATHMM-XF | HGMD | no | 92% | 0.914 [0.901–0.927] | -0.072 [-0.084, -0.061] |
| MutationAssessor | none | no | 92% | 0.912 [0.898–0.926] | -0.073 [-0.086, -0.061] |
| PROVEAN | none | no | 95% | 0.910 [0.896–0.924] | -0.076 [-0.088, -0.064] |
| popEVE | none (population) | no | 83% | 0.906 [0.892–0.921] | -0.079 [-0.091, -0.067] |
| SIFT4G | none | no | 95% | 0.898 [0.884–0.913] | -0.087 [-0.101, -0.074] |
| LIST-S2 | few parameters | no | 98% | 0.893 [0.880–0.906] | -0.092 [-0.105, -0.079] |
| PolyPhen-2 | UniProt | no | 93% | 0.884 [0.868–0.899] | -0.101 [-0.117, -0.087] |
| PrimateAI | none (population) | no | 98% | 0.877 [0.861–0.892] | -0.108 [-0.123, -0.095] |
| SIFT | none | no | 98% | 0.862 [0.845–0.878] | -0.123 [-0.138, -0.109] |
| phyloP (100 vertebrates) | none | no | 100% | 0.859 [0.843–0.874] | -0.126 [-0.142, -0.112] |
| DANN | none (simulated) | no | 100% | 0.827 [0.808–0.844] | -0.159 [-0.177, -0.142] |
| MPC | ClinVar (402 variants) | no | 86% | 0.817 [0.796–0.837] | -0.169 [-0.191, -0.149] |

## Rare variants

| Tool | Clinical training labels | Uses allele frequency | Coverage | AUC | Δ vs ours |
|---|---|---|---|---|---|
| MetaRNN | ClinVar | yes | 99% | 0.983 [0.975–0.991] | +0.009 [+0.001, +0.017] |
| **MissensExplorer v2.5 + frequency as a feature (exploratory)** | ClinVar (other genes) | yes | 100% | 0.983 [0.973–0.991] | +0.009 [+0.003, +0.016] |
| **MissensExplorer v2.5** | ClinVar (other genes) | no | 100% | 0.974 [0.961–0.985] | — |
| **MissensExplorer v2.5 + frequency rule (BA1)** | ClinVar (other genes) | yes | 100% | 0.974 [0.960–0.985] | +0.000 [+0.000, +0.000] |
| MutationTaster2021 | ClinVar, HGMD (+ look-up) | no | 96% | 0.974 [0.958–0.987] | -0.002 [-0.017, +0.012] |
| ClinPred | ClinVar | yes | 99% | 0.973 [0.958–0.985] | -0.002 [-0.011, +0.006] |
| BayesDel (with AF) | ClinVar, UniProt | yes | 100% | 0.972 [0.955–0.986] | -0.002 [-0.014, +0.007] |
| REVEL | HGMD | no | 96% | 0.972 [0.958–0.983] | -0.004 [-0.008, +0.001] |
| gMVP | ClinVar, HGMD, UniProt | no | 93% | 0.966 [0.949–0.979] | -0.008 [-0.018, +0.002] |
| BayesDel (no AF) | ClinVar, UniProt | no | 100% | 0.962 [0.944–0.976] | -0.013 [-0.025, -0.003] |
| VARITY_R | ClinVar, HGMD | no | 91% | 0.960 [0.945–0.973] | -0.013 [-0.022, -0.005] |
| VARITY_R (leave-one-out) | ClinVar, HGMD | no | 91% | 0.957 [0.941–0.970] | -0.017 [-0.026, -0.007] |
| PHACTboost | ClinVar | no | 93% | 0.957 [0.936–0.974] | -0.018 [-0.032, -0.005] |
| AlphaMissense | none (population) | no | 99% | 0.953 [0.938–0.967] | -0.022 [-0.033, -0.012] |
| MutPred2 | HGMD, SwissVar | no | 100% | 0.946 [0.928–0.963] | -0.028 [-0.042, -0.015] |
| MVP | ClinVar, HGMD, UniProt | no | 89% | 0.938 [0.915–0.958] | -0.037 [-0.057, -0.021] |
| DEOGEN2 | UniProt | no | 97% | 0.930 [0.906–0.949] | -0.044 [-0.059, -0.031] |
| MutationAssessor | none | no | 92% | 0.923 [0.904–0.941] | -0.050 [-0.065, -0.037] |
| ESM1b | none | no | 100% | 0.920 [0.898–0.941] | -0.054 [-0.074, -0.036] |
| Eigen | none | no | 92% | 0.916 [0.894–0.938] | -0.058 [-0.076, -0.041] |
| MetaSVM | UniProt | yes | 99% | 0.916 [0.887–0.940] | -0.058 [-0.083, -0.038] |
| gnomAD frequency alone (reference) | — | yes | 100% | 0.909 [0.890–0.927] | -0.065 [-0.086, -0.046] |
| MetaLR | UniProt | yes | 99% | 0.909 [0.882–0.931] | -0.066 [-0.087, -0.047] |
| PROVEAN | none | no | 95% | 0.905 [0.879–0.930] | -0.070 [-0.092, -0.050] |
| EVE | none | no | 59% | 0.901 [0.872–0.930] | -0.065 [-0.093, -0.041] |
| SIFT4G | none | no | 95% | 0.897 [0.873–0.920] | -0.077 [-0.099, -0.056] |
| FATHMM-XF | HGMD | no | 92% | 0.897 [0.868–0.922] | -0.077 [-0.101, -0.057] |
| MutFormer | HGMD | no | 99% | 0.889 [0.861–0.915] | -0.085 [-0.111, -0.060] |
| CADD | none (simulated) | no | 100% | 0.885 [0.857–0.913] | -0.089 [-0.114, -0.063] |
| PolyPhen-2 | UniProt | no | 93% | 0.883 [0.858–0.910] | -0.092 [-0.113, -0.070] |
| popEVE | none (population) | no | 83% | 0.856 [0.821–0.890] | -0.117 [-0.145, -0.087] |
| LIST-S2 | few parameters | no | 98% | 0.851 [0.822–0.880] | -0.123 [-0.151, -0.096] |
| phyloP (100 vertebrates) | none | no | 100% | 0.835 [0.803–0.864] | -0.139 [-0.167, -0.114] |
| SIFT | none | no | 98% | 0.805 [0.765–0.842] | -0.170 [-0.206, -0.136] |
| PrimateAI | none (population) | no | 98% | 0.793 [0.756–0.831] | -0.180 [-0.213, -0.146] |
| DANN | none (simulated) | no | 100% | 0.787 [0.746–0.825] | -0.187 [-0.224, -0.151] |
| MPC | ClinVar (402 variants) | no | 86% | 0.746 [0.709–0.786] | -0.229 [-0.263, -0.189] |

## New since 2023 variants

| Tool | Clinical training labels | Uses allele frequency | Coverage | AUC | Δ vs ours |
|---|---|---|---|---|---|
| MetaRNN | ClinVar | yes | 99% | 0.983 [0.972–0.992] | +0.012 [+0.003, +0.024] |
| **MissensExplorer v2.5 + frequency as a feature (exploratory)** | ClinVar (other genes) | yes | 100% | 0.982 [0.970–0.992] | +0.011 [+0.002, +0.021] |
| MutationTaster2021 | ClinVar, HGMD (+ look-up) | no | 96% | 0.972 [0.952–0.989] | +0.001 [-0.017, +0.018] |
| ClinPred | ClinVar | yes | 99% | 0.971 [0.953–0.986] | +0.001 [-0.011, +0.012] |
| **MissensExplorer v2.5 + frequency rule (BA1)** | ClinVar (other genes) | yes | 100% | 0.971 [0.953–0.985] | +0.000 [+0.000, +0.001] |
| **MissensExplorer v2.5** | ClinVar (other genes) | no | 100% | 0.971 [0.952–0.985] | — |
| BayesDel (with AF) | ClinVar, UniProt | yes | 100% | 0.969 [0.947–0.984] | -0.002 [-0.018, +0.010] |
| REVEL | HGMD | no | 96% | 0.966 [0.948–0.981] | -0.005 [-0.013, +0.003] |
| AlphaMissense | none (population) | no | 99% | 0.957 [0.938–0.973] | -0.014 [-0.027, -0.004] |
| PHACTboost | ClinVar | no | 93% | 0.956 [0.936–0.974] | -0.015 [-0.027, -0.004] |
| gMVP | ClinVar, HGMD, UniProt | no | 93% | 0.956 [0.934–0.975] | -0.017 [-0.033, -0.004] |
| VARITY_R | ClinVar, HGMD | no | 91% | 0.954 [0.933–0.972] | -0.016 [-0.032, -0.002] |
| VARITY_R (leave-one-out) | ClinVar, HGMD | no | 91% | 0.953 [0.933–0.972] | -0.017 [-0.031, -0.004] |
| BayesDel (no AF) | ClinVar, UniProt | no | 100% | 0.950 [0.927–0.969] | -0.021 [-0.037, -0.008] |
| MutPred2 | HGMD, SwissVar | no | 100% | 0.942 [0.919–0.963] | -0.029 [-0.047, -0.014] |
| gnomAD frequency alone (reference) | — | yes | 100% | 0.938 [0.916–0.957] | -0.033 [-0.055, -0.010] |
| MVP | ClinVar, HGMD, UniProt | no | 89% | 0.929 [0.895–0.956] | -0.042 [-0.069, -0.019] |
| MutationAssessor | none | no | 92% | 0.923 [0.896–0.950] | -0.046 [-0.068, -0.027] |
| Eigen | none | no | 92% | 0.921 [0.895–0.945] | -0.052 [-0.073, -0.033] |
| MutFormer | HGMD | no | 99% | 0.917 [0.888–0.943] | -0.054 [-0.079, -0.029] |
| DEOGEN2 | UniProt | no | 97% | 0.912 [0.881–0.940] | -0.058 [-0.081, -0.036] |
| FATHMM-XF | HGMD | no | 92% | 0.910 [0.879–0.936] | -0.063 [-0.091, -0.042] |
| ESM1b | none | no | 100% | 0.909 [0.878–0.937] | -0.062 [-0.089, -0.037] |
| MetaSVM | UniProt | yes | 99% | 0.907 [0.874–0.936] | -0.064 [-0.093, -0.038] |
| MetaLR | UniProt | yes | 99% | 0.906 [0.876–0.934] | -0.065 [-0.089, -0.041] |
| EVE | none | no | 59% | 0.906 [0.867–0.941] | -0.059 [-0.092, -0.029] |
| PolyPhen-2 | UniProt | no | 93% | 0.896 [0.863–0.926] | -0.077 [-0.107, -0.048] |
| PROVEAN | none | no | 95% | 0.893 [0.858–0.923] | -0.077 [-0.108, -0.051] |
| CADD | none (simulated) | no | 100% | 0.890 [0.857–0.921] | -0.081 [-0.108, -0.055] |
| SIFT4G | none | no | 95% | 0.886 [0.852–0.917] | -0.083 [-0.114, -0.053] |
| popEVE | none (population) | no | 83% | 0.867 [0.826–0.907] | -0.098 [-0.130, -0.066] |
| LIST-S2 | few parameters | no | 98% | 0.864 [0.828–0.896] | -0.107 [-0.139, -0.076] |
| phyloP (100 vertebrates) | none | no | 100% | 0.856 [0.820–0.890] | -0.115 [-0.146, -0.085] |
| SIFT | none | no | 98% | 0.812 [0.769–0.856] | -0.158 [-0.195, -0.119] |
| DANN | none (simulated) | no | 100% | 0.810 [0.763–0.854] | -0.161 [-0.203, -0.120] |
| PrimateAI | none (population) | no | 98% | 0.802 [0.754–0.841] | -0.169 [-0.208, -0.132] |
| MPC | ClinVar (402 variants) | no | 86% | 0.734 [0.682–0.785] | -0.236 [-0.286, -0.188] |

## Absent from gnomAD variants

| Tool | Clinical training labels | Uses allele frequency | Coverage | AUC | Δ vs ours |
|---|---|---|---|---|---|
| FATHMM-XF | HGMD | no | 92% | 1.000 [1.000–1.000] | +0.000 [+0.000, +0.000] |
| Eigen | none | no | 92% | 1.000 [1.000–1.000] | +0.000 [+0.000, +0.000] |
| DANN | none (simulated) | no | 100% | 0.980 [0.950–1.000] | +0.035 [+0.000, +0.089] |
| SIFT4G | none | no | 95% | 0.971 [0.939–0.996] | +0.027 [-0.013, +0.080] |
| AlphaMissense | none (population) | no | 99% | 0.961 [0.907–1.000] | +0.015 [-0.004, +0.045] |
| EVE | none | no | 59% | 0.956 [0.892–1.000] | -0.002 [-0.026, +0.015] |
| **MissensExplorer v2.5** | ClinVar (other genes) | no | 100% | 0.946 [0.876–1.000] | — |
| **MissensExplorer v2.5 + frequency rule (BA1)** | ClinVar (other genes) | yes | 100% | 0.946 [0.874–1.000] | +0.000 [+0.000, +0.000] |
| **MissensExplorer v2.5 + frequency as a feature (exploratory)** | ClinVar (other genes) | yes | 100% | 0.946 [0.874–1.000] | +0.000 [+0.000, +0.000] |
| MutationAssessor | none | no | 92% | 0.929 [0.839–0.998] | -0.016 [-0.051, +0.002] |
| CADD | none (simulated) | no | 100% | 0.923 [0.827–1.000] | -0.023 [-0.066, +0.000] |
| REVEL | HGMD | no | 96% | 0.875 [0.726–1.000] | -0.071 [-0.164, +0.000] |
| VARITY_R | ClinVar, HGMD | no | 91% | 0.867 [0.706–1.000] | -0.081 [-0.185, +0.000] |
| VARITY_R (leave-one-out) | ClinVar, HGMD | no | 91% | 0.864 [0.698–1.000] | -0.084 [-0.198, +0.000] |
| ESM1b | none | no | 100% | 0.836 [0.652–0.998] | -0.110 [-0.239, -0.002] |
| popEVE | none (population) | no | 83% | 0.835 [0.687–0.971] | -0.115 [-0.213, -0.029] |
| MetaRNN | ClinVar | yes | 99% | 0.808 [0.588–1.000] | -0.137 [-0.305, +0.000] |
| ClinPred | ClinVar | yes | 99% | 0.807 [0.584–1.000] | -0.139 [-0.302, +0.000] |
| LIST-S2 | few parameters | no | 98% | 0.802 [0.618–0.970] | -0.146 [-0.282, -0.030] |
| PHACTboost | ClinVar | no | 93% | 0.791 [0.551–1.000] | -0.161 [-0.351, +0.000] |
| MutPred2 | HGMD, SwissVar | no | 100% | 0.787 [0.544–1.000] | -0.159 [-0.348, +0.000] |
| MPC | ClinVar (402 variants) | no | 86% | 0.748 [0.489–0.987] | -0.198 [-0.405, -0.013] |
| PolyPhen-2 | UniProt | no | 93% | 0.729 [0.454–0.980] | -0.211 [-0.431, -0.020] |
| SIFT | none | no | 98% | 0.709 [0.408–1.000] | -0.237 [-0.490, +0.000] |
| gMVP | ClinVar, HGMD, UniProt | no | 93% | 0.672 [0.311–1.000] | -0.271 [-0.571, +0.000] |
| PrimateAI | none (population) | no | 98% | 0.669 [0.315–1.000] | -0.276 [-0.579, +0.000] |
| MetaLR | UniProt | yes | 99% | 0.656 [0.300–0.992] | -0.290 [-0.594, -0.008] |
| MutFormer | HGMD | no | 99% | 0.638 [0.596–0.680] | -0.308 [-0.397, -0.213] |
| phyloP (100 vertebrates) | none | no | 100% | 0.608 [0.203–0.992] | -0.338 [-0.696, -0.008] |
| DEOGEN2 | UniProt | no | 97% | 0.598 [0.184–0.996] | -0.347 [-0.711, -0.004] |
| MutationTaster2021 | ClinVar, HGMD (+ look-up) | no | 96% | 0.565 [0.522–0.612] | -0.329 [-0.369, -0.279] |
| MVP | ClinVar, HGMD, UniProt | no | 89% | 0.561 [0.106–1.000] | -0.385 [-0.792, +0.000] |
| PROVEAN | none | no | 95% | 0.530 [0.045–1.000] | -0.416 [-0.855, +0.000] |
| BayesDel (no AF) | ClinVar, UniProt | no | 100% | 0.525 [0.040–1.000] | -0.420 [-0.859, +0.000] |
| BayesDel (with AF) | ClinVar, UniProt | yes | 100% | 0.525 [0.042–1.000] | -0.420 [-0.859, +0.000] |
| gnomAD frequency alone (reference) | — | yes | 100% | 0.500 [0.500–0.500] | -0.446 [-0.500, -0.372] |
| MetaSVM | UniProt | yes | 99% | 0.495 [0.004–0.990] | -0.451 [-0.901, -0.010] |

## Grpmax below 0.01% variants

| Tool | Clinical training labels | Uses allele frequency | Coverage | AUC | Δ vs ours |
|---|---|---|---|---|---|
| MetaRNN | ClinVar | yes | 99% | 0.966 [0.942–0.987] | +0.032 [+0.007, +0.068] |
| MutationTaster2021 | ClinVar, HGMD (+ look-up) | no | 96% | 0.945 [0.895–0.983] | +0.007 [-0.035, +0.045] |
| REVEL | HGMD | no | 96% | 0.942 [0.898–0.977] | +0.005 [-0.008, +0.019] |
| **MissensExplorer v2.5 + frequency as a feature (exploratory)** | ClinVar (other genes) | yes | 100% | 0.940 [0.893–0.976] | +0.007 [-0.000, +0.015] |
| gMVP | ClinVar, HGMD, UniProt | no | 93% | 0.937 [0.888–0.972] | +0.005 [-0.028, +0.036] |
| BayesDel (with AF) | ClinVar, UniProt | yes | 100% | 0.936 [0.879–0.977] | +0.003 [-0.020, +0.022] |
| BayesDel (no AF) | ClinVar, UniProt | no | 100% | 0.933 [0.880–0.973] | +0.000 [-0.023, +0.021] |
| **MissensExplorer v2.5** | ClinVar (other genes) | no | 100% | 0.933 [0.879–0.972] | — |
| **MissensExplorer v2.5 + frequency rule (BA1)** | ClinVar (other genes) | yes | 100% | 0.933 [0.886–0.970] | +0.000 [+0.000, +0.000] |
| VARITY_R | ClinVar, HGMD | no | 91% | 0.930 [0.882–0.968] | +0.001 [-0.022, +0.026] |
| VARITY_R (leave-one-out) | ClinVar, HGMD | no | 91% | 0.926 [0.880–0.969] | -0.003 [-0.026, +0.021] |
| ClinPred | ClinVar | yes | 99% | 0.923 [0.872–0.963] | -0.010 [-0.034, +0.011] |
| MutPred2 | HGMD, SwissVar | no | 100% | 0.921 [0.867–0.963] | -0.012 [-0.053, +0.027] |
| PHACTboost | ClinVar | no | 93% | 0.916 [0.849–0.967] | -0.013 [-0.051, +0.020] |
| AlphaMissense | none (population) | no | 99% | 0.904 [0.852–0.944] | -0.030 [-0.058, -0.005] |
| MetaSVM | UniProt | yes | 99% | 0.893 [0.837–0.939] | -0.040 [-0.076, -0.002] |
| DEOGEN2 | UniProt | no | 97% | 0.890 [0.831–0.939] | -0.042 [-0.070, -0.015] |
| MVP | ClinVar, HGMD, UniProt | no | 89% | 0.888 [0.822–0.942] | -0.048 [-0.093, -0.008] |
| MutationAssessor | none | no | 92% | 0.888 [0.836–0.932] | -0.042 [-0.077, -0.011] |
| SIFT4G | none | no | 95% | 0.884 [0.838–0.924] | -0.054 [-0.100, -0.009] |
| Eigen | none | no | 92% | 0.880 [0.818–0.929] | -0.044 [-0.080, -0.009] |
| ESM1b | none | no | 100% | 0.875 [0.812–0.927] | -0.058 [-0.108, -0.014] |
| MetaLR | UniProt | yes | 99% | 0.871 [0.811–0.923] | -0.062 [-0.103, -0.022] |
| PROVEAN | none | no | 95% | 0.867 [0.795–0.923] | -0.068 [-0.123, -0.020] |
| PolyPhen-2 | UniProt | no | 93% | 0.865 [0.803–0.915] | -0.075 [-0.126, -0.029] |
| EVE | none | no | 59% | 0.840 [0.759–0.906] | -0.065 [-0.128, -0.001] |
| FATHMM-XF | HGMD | no | 92% | 0.839 [0.766–0.907] | -0.084 [-0.139, -0.035] |
| MutFormer | HGMD | no | 99% | 0.838 [0.772–0.905] | -0.095 [-0.160, -0.038] |
| CADD | none (simulated) | no | 100% | 0.834 [0.759–0.895] | -0.099 [-0.161, -0.048] |
| LIST-S2 | few parameters | no | 98% | 0.813 [0.746–0.871] | -0.120 [-0.187, -0.056] |
| DANN | none (simulated) | no | 100% | 0.772 [0.698–0.847] | -0.161 [-0.232, -0.085] |
| SIFT | none | no | 98% | 0.763 [0.687–0.842] | -0.175 [-0.245, -0.094] |
| phyloP (100 vertebrates) | none | no | 100% | 0.762 [0.684–0.827] | -0.171 [-0.227, -0.115] |
| popEVE | none (population) | no | 83% | 0.748 [0.635–0.837] | -0.183 [-0.267, -0.112] |
| gnomAD frequency alone (reference) | — | yes | 100% | 0.707 [0.638–0.772] | -0.226 [-0.303, -0.158] |
| PrimateAI | none (population) | no | 98% | 0.680 [0.594–0.767] | -0.251 [-0.330, -0.169] |
| MPC | ClinVar (402 variants) | no | 86% | 0.672 [0.580–0.766] | -0.264 [-0.353, -0.174] |

## Same variants for every tool (1290 variants scored by all tools with ≥ 90% coverage)

| Tool | All | Rare | New since 2023 |
|---|---|---|---|
| MetaRNN | 0.997 | 0.986 | 0.987 |
| MissensExplorer v2.5 + frequency as a feature (exploratory) | 0.996 | 0.982 | 0.984 |
| ClinPred | 0.994 | 0.973 | 0.978 |
| BayesDel (with AF) | 0.994 | 0.975 | 0.974 |
| MutationTaster2021 | 0.994 | 0.980 | 0.986 |
| MissensExplorer v2.5 + frequency rule (BA1) | 0.987 | 0.975 | 0.973 |
| MissensExplorer v2.5 | 0.987 | 0.975 | 0.973 |
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
| **MissensExplorer v2.5 (default verdict)** | 91.9% | 79.9% | 2.2% | 6.0% | 92.7% | 96.9% |
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
