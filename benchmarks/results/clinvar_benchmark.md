# ClinVar benchmark

1972 expert-reviewed (≥ 2 stars) germline missense variants from 1181 genes, GRCh38, sampled with seed 42 from ClinVar's variant_summary. ClinVar is not an input of the verdict. Run: MissensExplorer 2.3.0, Ensembl VEP REST, release 116, 2026-11-23T09:04:04+00:00.

The missense model (v2.3, trained 2026-11-23) was fitted on 9218 other ClinVar variants from 3263 genes; **no gene of this benchmark was used for training**, so these are held-out results (benchmarks/train_model.py). Its cut-offs (0.76/0.3) were chosen by cross-validation on the training genes, not on this benchmark.

Sensitivity = pathogenic variants called *damaging*; specificity = benign variants called *benign*; coverage = share of variants that get a decisive call (not *uncertain*/*not assessed*).

## All (n = 978 pathogenic, 994 benign)

| Method | Sensitivity | Specificity | Pathogenic called benign | Benign called damaging | Coverage | Accuracy when decisive |
|---|---|---|---|---|---|---|
| v2.3 default: model + frequency rule (BA1, gnomAD grpmax ≥ 5%) | 91.4% | 88.4% | 2.4% | 3.2% | 92.7% | 97.0% |
| v2.3 model alone (no frequency rule) | 91.4% | 88.1% | 2.4% | 3.2% | 92.5% | 97.0% |
| v2.3 model + frequency ≥ 1% | 91.4% | 91.4% | 2.4% | 2.4% | 93.8% | 97.5% |
| v2.2 model + frequency rule (BA1, 5%) | 93.1% | 88.4% | 2.4% | 5.1% | 94.5% | 96.0% |
| v2.1 vote | 89.1% | 87.3% | 2.7% | 5.2% | 92.1% | 95.7% |
| v2.1 vote + frequency rule (BA1, 5%) | 89.1% | 87.8% | 2.7% | 5.1% | 92.3% | 95.8% |
| ClinGen PP3/BP4 alone (REVEL, else AlphaMissense) | 88.2% | 76.7% | 1.4% | 3.3% | 84.8% | 97.2% |

AUC: v2.3 model score 0.985; v2.2 model score 0.985; v2.1 vote fraction 0.959 (coverage 99.4%); AlphaMissense 0.964 (coverage 98.7%); REVEL 0.982 (coverage 95.7%); CADD 0.924; EVE 0.928 (coverage 59.3%); popEVE 0.906 (coverage 82.8%); SIFT 0.862 (coverage 98.1%); PolyPhen-2 0.884 (coverage 93.1%).

| Predictor alone | Sensitivity | Specificity | Coverage | Accuracy when decisive |
|---|---|---|---|---|
| AlphaMissense | 82.4% | 88.2% | 92.1% | 92.7% |
| REVEL | 85.8% | 73.2% | 81.7% | 97.3% |
| EVE | 58.8% | 29.0% | 48.5% | 90.2% |
| CADD | 78.3% | 73.1% | 84.2% | 89.9% |
| PolyPhen | 67.7% | 63.6% | 74.9% | 87.6% |
| SIFT | 95.5% | 58.6% | 98.1% | 78.3% |

## Rare (gnomAD grpmax < 0.1% or absent) (n = 964 pathogenic, 199 benign)

| Method | Sensitivity | Specificity | Pathogenic called benign | Benign called damaging | Coverage | Accuracy when decisive |
|---|---|---|---|---|---|---|
| v2.3 default: model + frequency rule (BA1, gnomAD grpmax ≥ 5%) | 91.6% | 79.4% | 2.3% | 6.0% | 92.4% | 96.8% |
| v2.3 model alone (no frequency rule) | 91.6% | 79.4% | 2.3% | 6.0% | 92.4% | 96.8% |
| v2.3 model + frequency ≥ 1% | 91.6% | 79.4% | 2.3% | 6.0% | 92.4% | 96.8% |
| v2.2 model + frequency rule (BA1, 5%) | 93.4% | 80.4% | 2.3% | 10.6% | 94.8% | 96.1% |
| v2.1 vote | 89.2% | 81.9% | 2.6% | 9.0% | 91.7% | 96.0% |
| v2.1 vote + frequency rule (BA1, 5%) | 89.2% | 81.9% | 2.6% | 9.0% | 91.7% | 96.0% |
| ClinGen PP3/BP4 alone (REVEL, else AlphaMissense) | 88.5% | 71.4% | 1.5% | 6.0% | 87.8% | 97.5% |

AUC: v2.3 model score 0.974; v2.2 model score 0.973; v2.1 vote fraction 0.937; AlphaMissense 0.953 (coverage 99.3%); REVEL 0.972 (coverage 96.6%); CADD 0.885; EVE 0.901 (coverage 75.9%); popEVE 0.856 (coverage 86.8%); SIFT 0.805 (coverage 98.6%); PolyPhen-2 0.884 (coverage 91.0%).

## Evaluated since 2024 (n = 773 pathogenic, 769 benign)

| Method | Sensitivity | Specificity | Pathogenic called benign | Benign called damaging | Coverage | Accuracy when decisive |
|---|---|---|---|---|---|---|
| v2.3 default: model + frequency rule (BA1, gnomAD grpmax ≥ 5%) | 90.6% | 88.3% | 2.7% | 3.4% | 92.5% | 96.7% |
| v2.3 model alone (no frequency rule) | 90.6% | 88.0% | 2.7% | 3.4% | 92.3% | 96.7% |
| v2.3 model + frequency ≥ 1% | 90.6% | 90.6% | 2.7% | 2.7% | 93.3% | 97.1% |
| v2.2 model + frequency rule (BA1, 5%) | 92.5% | 88.3% | 2.7% | 5.3% | 94.4% | 95.7% |
| v2.1 vote | 88.9% | 87.6% | 2.7% | 5.1% | 92.2% | 95.8% |
| v2.1 vote + frequency rule (BA1, 5%) | 88.9% | 88.0% | 2.7% | 5.1% | 92.3% | 95.8% |
| ClinGen PP3/BP4 alone (REVEL, else AlphaMissense) | 87.7% | 75.4% | 1.7% | 3.4% | 84.1% | 97.0% |

AUC: v2.3 model score 0.984; v2.2 model score 0.983; v2.1 vote fraction 0.960; AlphaMissense 0.962 (coverage 98.9%); REVEL 0.980 (coverage 95.8%); CADD 0.921; EVE 0.931 (coverage 62.5%); popEVE 0.898 (coverage 81.8%); SIFT 0.862 (coverage 98.1%); PolyPhen-2 0.888 (coverage 92.9%).

## Model cut-offs (all variants, no frequency rule)

| Method | Sensitivity | Specificity | Pathogenic called benign | Benign called damaging | Coverage | Accuracy when decisive |
|---|---|---|---|---|---|---|
| 0.60/0.40 | 94.5% | 90.8% | 3.3% | 6.5% | 97.6% | 95.0% |
| 0.70/0.30 | 92.9% | 88.1% | 2.4% | 4.8% | 94.1% | 96.2% |
| 0.76/0.30 (default) | 91.4% | 88.1% | 2.4% | 3.2% | 92.5% | 97.0% |
| 0.80/0.20 | 89.5% | 82.8% | 1.2% | 2.6% | 88.0% | 97.8% |

## Model cut-offs (rare variants, no frequency rule)

| Method | Sensitivity | Specificity | Pathogenic called benign | Benign called damaging | Coverage | Accuracy when decisive |
|---|---|---|---|---|---|---|
| 0.60/0.40 | 94.7% | 83.9% | 3.1% | 13.1% | 97.7% | 95.1% |
| 0.70/0.30 | 93.2% | 79.4% | 2.3% | 10.1% | 94.4% | 96.2% |
| 0.76/0.30 (default) | 91.6% | 79.4% | 2.3% | 6.0% | 92.4% | 96.8% |
| 0.80/0.20 | 89.6% | 74.9% | 1.2% | 5.5% | 89.1% | 97.8% |


Caveats: REVEL was trained on disease-variant databases and AlphaMissense's thresholds were calibrated on ClinVar, so absolute figures are optimistic; the comparison between methods is the robust part. The frequency rule looks better here than it will in practice because ClinVar curators use the same gnomAD frequencies to classify variants as benign; it also mislabels pathogenic variants of recessive or founder conditions (that is why the default is 5% and not 1%).
