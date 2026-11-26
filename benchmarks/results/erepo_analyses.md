# eRepo comparison: secondary analyses

6,310 classifications, 48 panels, 133 genes. Produced by `benchmarks/erepo_analyses.py` from `erepo_validation_variants.csv.gz`.

## Agreement with panel-clustered 95% intervals

| Criterion | Positive agreement (panel-cluster CI) | Negative agreement (panel-cluster CI) |
|---|---|---|
| PP3 | 89.6% (86.8%–92.3%) | 78.0% (71.6%–82.2%) |
| BP4 | 76.7% (71.1%–87.5%) | 81.2% (69.4%–87.9%) |
| PM2 | 48.3% (40.0%–56.3%) | 99.2% (98.4%–99.9%) |
| PS1 | 60.0% (43.1%–74.1%) | 98.2% (96.6%–99.1%) |
| PM5 | 88.5% (83.6%–91.8%) | 79.1% (67.0%–87.7%) |
| BA1 | 21.2% (13.3%–30.1%) | 100.0% (100.0%–100.0%) |

## AUC against the panel's class: all

n = 3718 (2584 P/LP, 1134 B/LB), 125 genes.
- AUC: {"model_score": 0.961437, "revel": 0.95297, "alphamissense": 0.934466}; within genes (32,866 pairs): {"model_score": 0.959624, "revel": 0.942813, "alphamissense": 0.95255}
- model - revel: +0.0085 (DeLong +0.0040 to +0.0130, p = 0.00022; gene-cluster bootstrap -0.0011 to +0.0192)
- model - alphamissense: +0.0270 (DeLong +0.0218 to +0.0321, p = 1.5e-24; gene-cluster bootstrap +0.0156 to +0.0395)

## AUC against the panel's class: gene never seen

n = 3408 (2426 P/LP, 982 B/LB), 95 genes.
- AUC: {"model_score": 0.960812, "revel": 0.95177, "alphamissense": 0.936483}; within genes (31,343 pairs): {"model_score": 0.958492, "revel": 0.941087, "alphamissense": 0.951935}
- model - revel: +0.0090 (DeLong +0.0043 to +0.0138, p = 0.00018; gene-cluster bootstrap -0.0005 to +0.0194)
- model - alphamissense: +0.0243 (DeLong +0.0190 to +0.0297, p = 7.3e-19; gene-cluster bootstrap +0.0122 to +0.0361)

## PP3/BP4 strength

- Classifications with PP3 or BP4 evaluated: 5,656; exact agreement of signed points 23.6%; linear-weighted kappa 0.3137; kappa of direction 0.5611.
- Panel applied PP3: {"n": 3098, "same": 522, "higher": 2163, "lower_not_none": 92, "none": 321}; strengths {"1": 2502, "2": 510, "4": 86}; panels applying PP3 above supporting: 17.
- Panel applied BP4: {"n": 1011, "same": 174, "higher": 600, "lower_not_none": 1, "none": 236}.

## Points

- Panel class rebuilt from its met codes: 5,925 of 6,310 reproduce it.
- After replacing the panel's PP3/BP4 by the tool's (among those): {"Benign -> Likely Benign": 9, "Likely Benign -> Benign": 95, "Likely Benign -> Uncertain Significance": 13, "Likely Pathogenic -> Pathogenic": 446, "Likely Pathogenic -> Uncertain Significance": 59, "Pathogenic -> Likely Pathogenic": 38, "Uncertain Significance -> Likely Benign": 309, "Uncertain Significance -> Likely Pathogenic": 425}; by panel class: {"Uncertain Significance": 2210, "Likely Pathogenic": 1331, "Pathogenic": 1219, "Benign": 607, "Likely Benign": 558}.
- Class from the tool's automated codes alone (panel | automated): {"Benign | Benign": 97, "Benign | Likely Benign": 356, "Benign | Likely Pathogenic": 1, "Benign | Uncertain Significance": 160, "Likely Benign | Likely Benign": 381, "Likely Benign | Likely Pathogenic": 1, "Likely Benign | Uncertain Significance": 198, "Likely Pathogenic | Likely Benign": 18, "Likely Pathogenic | Likely Pathogenic": 326, "Likely Pathogenic | Uncertain Significance": 1029, "Pathogenic | Likely Benign": 9, "Pathogenic | Likely Pathogenic": 385, "Pathogenic | Uncertain Significance": 870, "Uncertain Significance | Likely Benign": 619, "Uncertain Significance | Likely Pathogenic": 96, "Uncertain Significance | Uncertain Significance": 1764}.

## 'Not met' lists

- Classifications listing >= 15 codes as not met: 1,484; panels where they are the majority: {"Brain Malformations VCEP": 0.609, "CDH1 VCEP": 0.889, "DICER1 and miRNA-Processing Gene VCEP": 0.576, "FBN1 VCEP": 0.517, "Familial Hypercholesterolemia VCEP": 0.633, "Myeloid Malignancy VCEP": 0.973}.
- Agreement without them: PP3 89.2% / 75.9%; BP4 82.2% / 74.3%; PM2 45.5% / 99.6%; PS1 56.9% / 96.8%; PM5 89.2% / 70.3%; BA1 22.2% / 100.0%
- BP4 specificity without classifications where the panel applied BA1 or BS1: 87.8%.

## PM2 and PS1/PM5

- PM2 missed: 2,146; by grpmax {"no frequency": 1, "< 1e-5": 810, "1e-5 to 1e-4": 1007, ">= 1e-4": 328}; median grpmax 2.236e-05; sensitivity by period {"2021 or earlier": 0.405458, "2022-2023": 0.448737, "2024 or later": 0.503714}.
- PS1/PM5 applied by the tool, by the reference's stars and class: {"PS1": {"1 stars, P": {"panel met": 1, "panel not met": 6}, "2 stars, LP": {"panel met": 0, "panel not met": 1}, "2 stars, P": {"panel met": 3, "panel not met": 7}, "3 stars, LP": {"panel met": 12, "panel not met": 11}, "3 stars, P": {"panel met": 23, "panel not met": 10}}, "PM5": {"1 stars, P": {"panel met": 23, "panel not met": 50}, "2 stars, LP": {"panel met": 2, "panel not met": 27}, "2 stars, P": {"panel met": 66, "panel not met": 135}, "3 stars, LP": {"panel met": 159, "panel not met": 189}, "3 stars, P": {"panel met": 391, "panel not met": 126}}}.

## Incorporation bias

Variants whose rebuilt class does not depend on PP3/BP4: 2799 (1725 P/LP, 1074 B/LB). AUC {"model_score": 0.956672, "revel": 0.945324, "alphamissense": 0.928889}; model - REVEL +0.0113 (DeLong +0.0062 to +0.0165; gene-cluster -0.0005 to +0.0245).

BP4 given with SpliceAI above 0.1: 0.
