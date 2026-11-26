# PP3/BP4 bands by missing predictor scores

Produced by `benchmarks/calibration_by_inputs.py`.

Bands whose point estimate / lower Katz bound reach the required ratio, by test set:

- ClinVar 2020 set: 8 / 6 of 8
- ClinVar 2020 set, shipped model, training variants removed: 8 / 5 of 8
- held-out benchmark, rare variants: 8 / 3 of 8
- held-out benchmark, rare variants labelled since 2023: 7 / 1 of 8

## ClinVar 2020 set, shipped model: all four scores (2068 P / 4002 B)

| Band | P | B | LR (Katz 95%) | Required |
|---|---|---|---|---|
| 4 | 1062 | 47 | 43.73 (32.81–58.28) | 33.53 |
| 3 | 182 | 26 | 13.55 (9.01–20.36) | 13.93 |
| 2 | 264 | 92 | 5.55 (4.41–7.0) | 5.79 |
| 1 | 231 | 138 | 3.24 (2.64–3.97) | 2.41 |
| -1 | 54 | 497 | 4.76 (3.61–6.27) | 2.41 |
| -2 | 54 | 1251 | 11.97 (9.16–15.64) | 5.79 |
| -3 | 14 | 1005 | 37.09 (21.95–62.69) | 13.93 |
| -4 | 1 | 249 | 128.67 (18.07–916.4) | 33.53 |

## ClinVar 2020 set, shipped model: popEVE missing (AlphaMissense and REVEL present) (378 P / 1380 B)

| Band | P | B | LR (Katz 95%) | Required |
|---|---|---|---|---|
| 4 | 124 | 6 | 75.45 (33.52–169.83) | 33.53 |
| 3 | 34 | 3 | 41.38 (12.78–133.98) | 13.93 |
| 2 | 66 | 21 | 11.47 (7.12–18.5) | 5.79 |
| 1 | 53 | 44 | 4.4 (3.0–6.45) | 2.41 |
| -1 | 16 | 214 | 3.66 (2.23–6.01) | 2.41 |
| -2 | 12 | 456 | 10.41 (5.93–18.25) | 5.79 |
| -3 | 2 | 363 | 49.72 (12.44–198.62) | 13.93 |
| -4 | 0 | 30 | ∞ (1.03–273.15) | 33.53 |

## ClinVar 2020 set, shipped model: AlphaMissense or REVEL missing (115 P / 390 B)

| Band | P | B | LR (Katz 95%) | Required |
|---|---|---|---|---|
| 4 | 17 | 0 | ∞ (7.15–1946.9) | 33.53 |
| 3 | 4 | 0 | ∞ (1.65–559.35) | 13.93 |
| 2 | 25 | 0 | ∞ (10.55–2802.21) | 5.79 |
| 1 | 28 | 8 | 11.87 (5.56–25.32) | 2.41 |
| -1 | 4 | 130 | 9.58 (3.62–25.36) | 2.41 |
| -2 | 3 | 116 | 11.4 (3.69–35.19) | 5.79 |
| -3 | 0 | 40 | ∞ (1.49–387.82) | 13.93 |
| -4 | 0 | 1 | ∞ (0.04–21.7) | 33.53 |
