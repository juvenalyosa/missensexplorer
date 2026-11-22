# ClinVar benchmark

1972 expert-reviewed (≥ 2 stars) germline missense variants from 1181 genes, GRCh38, sampled with seed 42 from ClinVar's variant_summary. ClinVar is not used by the verdict. Run: MissensExplorer 2.1.0, Ensembl VEP REST, release 116, 2026-11-22T17:59:16+00:00.

Sensitivity = pathogenic variants called *damaging*; specificity = benign variants called *benign*; coverage = share of variants that get a decisive call (not *uncertain*/*not assessed*).

## Verdict — all (n = 978 pathogenic, 994 benign)

| Sensitivity | Specificity | Pathogenic called benign | Benign called damaging | Coverage | Accuracy when decisive |
|---|---|---|---|---|---|
| 89.1% | 87.3% | 2.7% | 5.2% | 92.1% | 95.7% |

| Predictor | Sensitivity | Specificity | Coverage | Accuracy when decisive | AUC (score) |
|---|---|---|---|---|---|
| AlphaMissense | 82.4% | 88.2% | 92.1% | 92.7% | 0.964 |
| REVEL | 85.8% | 73.2% | 81.7% | 97.3% | 0.982 |
| EVE | 58.8% | 29.0% | 48.5% | 90.2% | 0.928 |
| CADD | 78.3% | 73.1% | 84.2% | 89.9% | 0.924 |
| PolyPhen-2 | 67.7% | 63.6% | 74.9% | 87.6% | 0.884 |
| SIFT | 95.5% | 58.6% | 98.1% | 78.3% | 0.862 |

AUC of the verdict's effect score: 0.959; popEVE: 0.906 (coverage 82.8%).

## Verdict — evaluated since 2024 (n = 773 pathogenic, 769 benign)

| Sensitivity | Specificity | Pathogenic called benign | Benign called damaging | Coverage | Accuracy when decisive |
|---|---|---|---|---|---|
| 88.9% | 87.6% | 2.7% | 5.1% | 92.2% | 95.8% |

| Predictor | Sensitivity | Specificity | Coverage | Accuracy when decisive | AUC (score) |
|---|---|---|---|---|---|
| AlphaMissense | 80.7% | 89.7% | 92.3% | 92.3% | 0.962 |
| REVEL | 85.8% | 71.8% | 81.2% | 97.0% | 0.980 |
| EVE | 59.0% | 34.6% | 51.6% | 90.8% | 0.931 |
| CADD | 77.5% | 73.0% | 83.9% | 89.7% | 0.921 |
| PolyPhen-2 | 69.2% | 63.2% | 75.1% | 88.2% | 0.888 |
| SIFT | 95.3% | 57.6% | 98.1% | 78.0% | 0.862 |

AUC of the verdict's effect score: 0.960; popEVE: 0.898 (coverage 81.8%).

## Alternative weights and cut-offs (2-fold cross-validation, split by gene)

| Scheme | Sensitivity | Specificity | Pathogenic called benign | Benign called damaging | Coverage | Accuracy when decisive |
|---|---|---|---|---|---|---|
| v2.0 (SIFT/PolyPhen always vote), cut-offs 0.70/0.30 | 90.8% | 81.8% | 1.9% | 6.1% | 90.3% | 95.5% |
| v2.0 (SIFT/PolyPhen always vote), cut-offs 0.70/0.30, min. weight 1.5 | 90.5% | 81.4% | 1.8% | 5.1% | 89.4% | 96.1% |
| v2.0 (SIFT/PolyPhen always vote), cut-offs 0.60/0.40 | 93.8% | 86.8% | 3.0% | 9.2% | 96.4% | 93.7% |
| v2.0 (SIFT/PolyPhen always vote), cut-offs 0.60/0.40, min. weight 1.5 | 93.5% | 86.3% | 2.9% | 8.2% | 95.3% | 94.2% |
| v2.0 (SIFT/PolyPhen always vote), cut-offs 0.80/0.20 | 84.4% | 77.9% | 1.5% | 5.3% | 84.5% | 95.9% |
| v2.0 (SIFT/PolyPhen always vote), cut-offs 0.80/0.20, min. weight 1.5 | 84.1% | 77.5% | 1.4% | 4.3% | 83.6% | 96.5% |
| equal weights, cut-offs 0.70/0.30 | 91.6% | 76.3% | 1.7% | 8.3% | 88.8% | 94.4% |
| equal weights, cut-offs 0.70/0.30, min. weight 1.5 | 91.5% | 76.1% | 1.7% | 8.0% | 88.5% | 94.5% |
| equal weights, cut-offs 0.60/0.40 | 95.0% | 83.3% | 3.3% | 12.2% | 96.9% | 91.9% |
| equal weights, cut-offs 0.60/0.40, min. weight 1.5 | 94.9% | 83.1% | 3.3% | 11.9% | 96.6% | 92.1% |
| equal weights, cut-offs 0.80/0.20 | 88.2% | 67.1% | 1.3% | 5.7% | 81.1% | 95.6% |
| equal weights, cut-offs 0.80/0.20, min. weight 1.5 | 88.1% | 66.9% | 1.3% | 5.4% | 80.8% | 95.8% |
| without SIFT/PolyPhen, cut-offs 0.70/0.30 | 89.1% | 87.3% | 2.7% | 5.2% | 92.1% | 95.7% |
| without SIFT/PolyPhen, cut-offs 0.70/0.30, min. weight 1.5 | 86.7% | 83.3% | 2.4% | 3.5% | 87.9% | 96.7% |
| without SIFT/PolyPhen, cut-offs 0.60/0.40 | 92.5% | 88.5% | 3.2% | 6.3% | 95.2% | 95.0% |
| without SIFT/PolyPhen, cut-offs 0.60/0.40, min. weight 1.5 | 90.1% | 84.5% | 3.0% | 4.5% | 91.0% | 95.9% |
| without SIFT/PolyPhen, cut-offs 0.80/0.20 | 84.2% | 84.3% | 1.8% | 4.9% | 87.6% | 96.1% |
| without SIFT/PolyPhen, cut-offs 0.80/0.20, min. weight 1.5 | 81.8% | 80.3% | 1.6% | 3.2% | 83.4% | 97.1% |
| SIFT/PolyPhen only as fallback (v2.1 default), cut-offs 0.70/0.30 | 89.1% | 87.3% | 2.7% | 5.2% | 92.1% | 95.7% |
| SIFT/PolyPhen only as fallback (v2.1 default), cut-offs 0.70/0.30, min. weight 1.5 | 86.7% | 83.3% | 2.4% | 3.5% | 87.9% | 96.7% |
| SIFT/PolyPhen only as fallback (v2.1 default), cut-offs 0.60/0.40 | 92.5% | 88.5% | 3.2% | 6.3% | 95.2% | 95.0% |
| SIFT/PolyPhen only as fallback (v2.1 default), cut-offs 0.60/0.40, min. weight 1.5 | 90.1% | 84.5% | 3.0% | 4.5% | 91.0% | 95.9% |
| SIFT/PolyPhen only as fallback (v2.1 default), cut-offs 0.80/0.20 | 84.2% | 84.3% | 1.8% | 4.9% | 87.6% | 96.1% |
| SIFT/PolyPhen only as fallback (v2.1 default), cut-offs 0.80/0.20, min. weight 1.5 | 81.8% | 80.3% | 1.6% | 3.2% | 83.4% | 97.1% |

Caveat: REVEL was trained on disease-variant databases and AlphaMissense's thresholds were calibrated on ClinVar, so these figures are optimistic; the *evaluated since 2024* rows are a less circular check.
