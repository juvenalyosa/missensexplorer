# The ten largest panels' own rules against their decisions (a consistency check)

Produced by `benchmarks/erepo_vcep_rules.py` (rules in `src/missensexplorer/vcep_rules.py`, encoded by hand from ClinGen CSpec; see the script for the limits). Intervals: bootstrap over the 10 panels.

| Criterion | Evaluated | Tool: positive / negative agreement | Panel rules: positive / negative agreement |
|---|---|---|---|
| PP3 | 2,824 | 91.8% (88.3–95.7) / 77.9% (64.5–84.9) | 98.4% (97.6–98.9) / 95.6% (81.0–99.7) |
| BP4 | 1,668 | 75.4% (70.0–91.7) / 86.7% (68.7–94.0) | 97.6% (86.8–100.0) / 97.1% (94.0–99.5) |
| PM2 | 2,878 | 47.4% (33.9–58.8) / 98.6% (97.7–100.0) | 89.1% (79.8–95.6) / 66.0% (43.1–96.6) |
| BA1 | 1,511 | 16.3% (3.5–31.9) / 100.0% (100.0–100.0) | 91.8% (83.2–100.0) / 98.1% (91.1–99.5) |

Exact agreement of the PP3/BP4 strength with the panels':

| Classifications | n | Tool | Tool capped at the specification's maximum | Panel rules |
|---|---|---|---|---|
| panel_pp3_met | 1,652 | 16.1% (10.5–21.8) | 79.9% (58.9–94.0) | 85.8% (61.4–98.5) |
| panel_pp3_or_bp4_met | 2,193 | 17.6% | 78.7% | 88.7% |
| pp3_or_bp4_evaluated | 3,013 | 25.2% (19.3–29.7) | 70.4% (59.0–80.3) | 88.9% (70.0–98.0) |

PM2 without the Myeloid Malignancy panel: panel rules 86.4% / 93.6%; tool 41.9% / 99.7%. Rule PM2 against the panel: 195 (174 Myeloid Malignancy; 176 in classifications listing 15 or more criteria as not met).

By approval period:

```
{
 "2016-2021": {
  "classifications": 559,
  "PP3": {
   "tool": {
    "evaluated": 449,
    "tp": 312,
    "fn": 25,
    "fp": 23,
    "tn": 89,
    "positive_agreement": 0.925816,
    "negative_agreement": 0.794643
   },
   "panel_rules": {
    "evaluated": 449,
    "tp": 331,
    "fn": 6,
    "fp": 40,
    "tn": 72,
    "positive_agreement": 0.982196,
    "negative_agreement": 0.642857
   }
  },
  "BP4": {
   "tool": {
    "evaluated": 130,
    "tp": 27,
    "fn": 2,
    "fp": 12,
    "tn": 89,
    "positive_agreement": 0.931034,
    "negative_agreement": 0.881188
   },
   "panel_rules": {
    "evaluated": 130,
    "tp": 22,
    "fn": 7,
    "fp": 6,
    "tn": 95,
    "positive_agreement": 0.758621,
    "negative_agreement": 0.940594
   }
  },
  "strength": {
   "n": 469,
   "tool": 0.238806,
   "tool_capped": 0.547974,
   "rule": 0.490405
  },
  "strength_panel_pp3_met": {
   "n": 337,
   "tool": 0.148368,
   "tool_capped": 0.522255,
   "rule": 0.448071
  }
 },
 "2022-2023": {
  "classifications": 730,
  "PP3": {
   "tool": {
    "evaluated": 644,
    "tp": 404,
    "fn": 35,
    "fp": 56,
    "tn": 149,
    "positive_agreement": 0.920273,
    "negative_agreement": 0.726829
   },
   "panel_rules": {
    "evaluated": 644,
    "tp": 436,
    "fn": 3,
    "fp": 7,
    "tn": 198,
    "positive_agreement": 0.993166,
    "negative_agreement": 0.965854
   }
  },
  "BP4": {
   "tool": {
    "evaluated": 390,
    "tp": 56,
    "fn": 10,
    "fp": 36,
    "tn": 288,
    "positive_agreement": 0.848485,
    "negative_agreement": 0.888889
   },
   "panel_rules": {
    "evaluated": 390,
    "tp": 65,
    "fn": 1,
    "fp": 12,
    "tn": 312,
    "positive_agreement": 0.984848,
    "negative_agreement": 0.962963
   }
  },
  "strength": {
   "n": 696,
   "tool": 0.275862,
   "tool_capped": 0.758621,
   "rule": 0.931034
  },
  "strength_panel_pp3_met": {
   "n": 439,
   "tool": 0.202733,
   "tool_capped": 0.865604,
   "rule": 0.940774
  }
 },
 "2024-2026": {
  "classifications": 2171,
  "PP3": {
   "tool": {
    "evaluated": 1731,
    "tp": 800,
    "fn": 76,
    "fp": 180,
    "tn": 675,
    "positive_agreement": 0.913242,
    "negative_agreement": 0.789474
   },
   "panel_rules": {
    "evaluated": 1731,
    "tp": 858,
    "fn": 18,
    "fp": 4,
    "tn": 851,
    "positive_agreement": 0.979452,
    "negative_agreement": 0.995322
   }
  },
  "BP4": {
   "tool": {
    "evaluated": 1148,
    "tp": 325,
    "fn": 121,
    "fp": 102,
    "tn": 600,
    "positive_agreement": 0.7287,
    "negative_agreement": 0.854701
   },
   "panel_rules": {
    "evaluated": 1148,
    "tp": 441,
    "fn": 5,
    "fp": 15,
    "tn": 687,
    "positive_agreement": 0.988789,
    "negative_agreement": 0.978632
   }
  },
  "strength": {
   "n": 1848,
   "tool": 0.245671,
   "tool_capped": 0.722944,
   "rule": 0.975108
  },
  "strength_panel_pp3_met": {
   "n": 876,
   "tool": 0.144977,
   "tool_capped": 0.872146,
   "rule": 0.974886
  }
 }
}
```

By panel:

```
{
 "Familial Hypercholesterolemia VCEP": {
  "classifications": 431,
  "PP3": {
   "tool": {
    "evaluated": 414,
    "tp": 219,
    "fn": 37,
    "fp": 11,
    "tn": 147,
    "positive_agreement": 0.855469,
    "negative_agreement": 0.93038
   },
   "panel_rules": {
    "evaluated": 414,
    "tp": 254,
    "fn": 2,
    "fp": 0,
    "tn": 158,
    "positive_agreement": 0.992188,
    "negative_agreement": 1.0
   }
  },
  "BP4": {
   "tool": {
    "evaluated": 346,
    "tp": 55,
    "fn": 3,
    "fp": 25,
    "tn": 263,
    "positive_agreement": 0.948276,
    "negative_agreement": 0.913194
   },
   "panel_rules": {
    "evaluated": 346,
    "tp": 58,
    "fn": 0,
    "fp": 2,
    "tn": 286,
    "positive_agreement": 1.0,
    "negative_agreement": 0.993056
   }
  },
  "PM2": {
   "tool": {
    "evaluated": 431,
    "tp": 171,
    "fn": 230,
    "fp": 0,
    "tn": 30,
    "positive_agreement": 0.426434,
    "negative_agreement": 1.0
   },
   "panel_rules": {
    "evaluated": 431,
    "tp": 386,
    "fn": 15,
    "fp": 4,
    "tn": 26,
    "positive_agreement": 0.962594,
    "negative_agreement": 0.866667
   }
  },
  "BA1": {
   "tool": {
    "evaluated": 306,
    "tp": 1,
    "fn": 3,
    "fp": 0,
    "tn": 302,
    "positive_agreement": 0.25,
    "negative_agreement": 1.0
   },
   "panel_rules": {
    "evaluated": 306,
    "tp": 4,
    "fn": 0,
    "fp": 0,
    "tn": 302,
    "positive_agreement": 1.0,
    "negative_agreement": 1.0
   }
  },
  "strength": {
   "n": 427,
   "tool": 0.316159,
   "tool_capped": 0.814988,
   "rule": 0.990632
  }
 },
 "Limb Girdle Muscular Dystrophy VCEP": {
  "classifications": 195,
  "PP3": {
   "tool": {
    "evaluated": 179,
    "tp": 130,
    "fn": 11,
    "fp": 13,
    "tn": 25,
    "positive_agreement": 0.921986,
    "negative_agreement": 0.657895
   },
   "panel_rules": {
    "evaluated": 179,
    "tp": 138,
    "fn": 3,
    "fp": 2,
    "tn": 36,
    "positive_agreement": 0.978723,
    "negative_agreement": 0.947368
   }
  },
  "BP4": {
   "tool": {
    "evaluated": 48,
    "tp": 1,
    "fn": 0,
    "fp": 15,
    "tn": 32,
    "positive_agreement": 1.0,
    "negative_agreement": 0.680851
   },
   "panel_rules": {
    "evaluated": 48,
    "tp": 1,
    "fn": 0,
    "fp": 0,
    "tn": 47,
    "positive_agreement": 1.0,
    "negative_agreement": 1.0
   }
  },
  "PM2": {
   "tool": {
    "evaluated": 174,
    "tp": 24,
    "fn": 88,
    "fp": 0,
    "tn": 62,
    "positive_agreement": 0.214286,
    "negative_agreement": 1.0
   },
   "panel_rules": {
    "evaluated": 174,
    "tp": 104,
    "fn": 8,
    "fp": 7,
    "tn": 55,
    "positive_agreement": 0.928571,
    "negative_agreement": 0.887097
   }
  },
  "BA1": {
   "tool": {
    "evaluated": 42,
    "tp": 2,
    "fn": 18,
    "fp": 0,
    "tn": 22,
    "positive_agreement": 0.1,
    "negative_agreement": 1.0
   },
   "panel_rules": {
    "evaluated": 42,
    "tp": 20,
    "fn": 0,
    "fp": 1,
    "tn": 21,
    "positive_agreement": 1.0,
    "negative_agreement": 0.954545
   }
  },
  "strength": {
   "n": 192,
   "tool": 0.239583,
   "tool_capped": 0.791667,
   "rule": 0.973958
  }
 },
 "Lysosomal Diseases VCEP": {
  "classifications": 232,
  "PP3": {
   "tool": {
    "evaluated": 204,
    "tp": 145,
    "fn": 30,
    "fp": 5,
    "tn": 24,
    "positive_agreement": 0.828571,
    "negative_agreement": 0.827586
   },
   "panel_rules": {
    "evaluated": 204,
    "tp": 172,
    "fn": 3,
    "fp": 0,
    "tn": 29,
    "positive_agreement": 0.982857,
    "negative_agreement": 1.0
   }
  },
  "BP4": {
   "tool": {
    "evaluated": 43,
    "tp": 16,
    "fn": 2,
    "fp": 2,
    "tn": 23,
    "positive_agreement": 0.888889,
    "negative_agreement": 0.92
   },
   "panel_rules": {
    "evaluated": 43,
    "tp": 18,
    "fn": 0,
    "fp": 1,
    "tn": 24,
    "positive_agreement": 1.0,
    "negative_agreement": 0.96
   }
  },
  "PM2": {
   "tool": {
    "evaluated": 214,
    "tp": 30,
    "fn": 175,
    "fp": 0,
    "tn": 9,
    "positive_agreement": 0.146341,
    "negative_agreement": 1.0
   },
   "panel_rules": {
    "evaluated": 214,
    "tp": 204,
    "fn": 1,
    "fp": 2,
    "tn": 7,
    "positive_agreement": 0.995122,
    "negative_agreement": 0.777778
   }
  },
  "BA1": {
   "tool": {
    "evaluated": 16,
    "tp": 8,
    "fn": 5,
    "fp": 0,
    "tn": 3,
    "positive_agreement": 0.615385,
    "negative_agreement": 1.0
   },
   "panel_rules": {
    "evaluated": 16,
    "tp": 13,
    "fn": 0,
    "fp": 0,
    "tn": 3,
    "positive_agreement": 1.0,
    "negative_agreement": 1.0
   }
  },
  "strength": {
   "n": 219,
   "tool": 0.273973,
   "tool_capped": 0.753425,
   "rule": 0.972603
  }
 },
 "Malignant Hyperthermia Susceptibility VCEP": {
  "classifications": 244,
  "PP3": {
   "tool": {
    "evaluated": 221,
    "tp": 128,
    "fn": 1,
    "fp": 49,
    "tn": 43,
    "positive_agreement": 0.992248,
    "negative_agreement": 0.467391
   },
   "panel_rules": {
    "evaluated": 221,
    "tp": 127,
    "fn": 2,
    "fp": 0,
    "tn": 92,
    "positive_agreement": 0.984496,
    "negative_agreement": 1.0
   }
  },
  "BP4": {
   "tool": {
    "evaluated": 138,
    "tp": 7,
    "fn": 4,
    "fp": 1,
    "tn": 126,
    "positive_agreement": 0.636364,
    "negative_agreement": 0.992126
   },
   "panel_rules": {
    "evaluated": 138,
    "tp": 11,
    "fn": 0,
    "fp": 0,
    "tn": 127,
    "positive_agreement": 1.0,
    "negative_agreement": 1.0
   }
  },
  "PM2": {
   "tool": {
    "evaluated": 0,
    "tp": 0,
    "fn": 0,
    "fp": 0,
    "tn": 0,
    "positive_agreement": null,
    "negative_agreement": null
   },
   "panel_rules": {
    "evaluated": 0,
    "tp": 0,
    "fn": 0,
    "fp": 0,
    "tn": 0,
    "positive_agreement": null,
    "negative_agreement": null
   }
  },
  "BA1": {
   "tool": {
    "evaluated": 127,
    "tp": 3,
    "fn": 12,
    "fp": 0,
    "tn": 112,
    "positive_agreement": 0.2,
    "negative_agreement": 1.0
   },
   "panel_rules": {
    "evaluated": 127,
    "tp": 15,
    "fn": 0,
    "fp": 0,
    "tn": 112,
    "positive_agreement": 1.0,
    "negative_agreement": 1.0
   }
  },
  "strength": {
   "n": 227,
   "tool": 0.378855,
   "tool_capped": 0.718062,
   "rule": 0.986784
  }
 },
 "Monogenic Diabetes VCEP": {
  "classifications": 437,
  "PP3": {
   "tool": {
    "evaluated": 418,
    "tp": 328,
    "fn": 21,
    "fp": 3,
    "tn": 66,
    "positive_agreement": 0.939828,
    "negative_agreement": 0.956522
   },
   "panel_rules": {
    "evaluated": 418,
    "tp": 344,
    "fn": 5,
    "fp": 0,
    "tn": 69,
    "positive_agreement": 0.985673,
    "negative_agreement": 1.0
   }
  },
  "BP4": {
   "tool": {
    "evaluated": 67,
    "tp": 4,
    "fn": 1,
    "fp": 19,
    "tn": 43,
    "positive_agreement": 0.8,
    "negative_agreement": 0.693548
   },
   "panel_rules": {
    "evaluated": 67,
    "tp": 5,
    "fn": 0,
    "fp": 0,
    "tn": 62,
    "positive_agreement": 1.0,
    "negative_agreement": 1.0
   }
  },
  "PM2": {
   "tool": {
    "evaluated": 408,
    "tp": 190,
    "fn": 169,
    "fp": 0,
    "tn": 49,
    "positive_agreement": 0.529248,
    "negative_agreement": 1.0
   },
   "panel_rules": {
    "evaluated": 408,
    "tp": 270,
    "fn": 89,
    "fp": 0,
    "tn": 49,
    "positive_agreement": 0.752089,
    "negative_agreement": 1.0
   }
  },
  "BA1": {
   "tool": {
    "evaluated": 39,
    "tp": 1,
    "fn": 18,
    "fp": 0,
    "tn": 20,
    "positive_agreement": 0.052632,
    "negative_agreement": 1.0
   },
   "panel_rules": {
    "evaluated": 39,
    "tp": 19,
    "fn": 0,
    "fp": 3,
    "tn": 17,
    "positive_agreement": 1.0,
    "negative_agreement": 0.85
   }
  },
  "strength": {
   "n": 423,
   "tool": 0.148936,
   "tool_capped": 0.86052,
   "rule": 0.98818
  }
 },
 "Myeloid Malignancy VCEP": {
  "classifications": 738,
  "PP3": {
   "tool": {
    "evaluated": 719,
    "tp": 146,
    "fn": 4,
    "fp": 112,
    "tn": 457,
    "positive_agreement": 0.973333,
    "negative_agreement": 0.803163
   },
   "panel_rules": {
    "evaluated": 719,
    "tp": 145,
    "fn": 5,
    "fp": 2,
    "tn": 567,
    "positive_agreement": 0.966667,
    "negative_agreement": 0.996485
   }
  },
  "BP4": {
   "tool": {
    "evaluated": 725,
    "tp": 235,
    "fn": 108,
    "fp": 17,
    "tn": 365,
    "positive_agreement": 0.685131,
    "negative_agreement": 0.955497
   },
   "panel_rules": {
    "evaluated": 725,
    "tp": 343,
    "fn": 0,
    "fp": 15,
    "tn": 367,
    "positive_agreement": 1.0,
    "negative_agreement": 0.960733
   }
  },
  "PM2": {
   "tool": {
    "evaluated": 736,
    "tp": 333,
    "fn": 159,
    "fp": 7,
    "tn": 237,
    "positive_agreement": 0.676829,
    "negative_agreement": 0.971311
   },
   "panel_rules": {
    "evaluated": 736,
    "tp": 487,
    "fn": 5,
    "fp": 174,
    "tn": 70,
    "positive_agreement": 0.989837,
    "negative_agreement": 0.286885
   }
  },
  "BA1": {
   "tool": {
    "evaluated": 718,
    "tp": 0,
    "fn": 12,
    "fp": 0,
    "tn": 706,
    "positive_agreement": 0.0,
    "negative_agreement": 1.0
   },
   "panel_rules": {
    "evaluated": 718,
    "tp": 7,
    "fn": 5,
    "fp": 5,
    "tn": 701,
    "positive_agreement": 0.583333,
    "negative_agreement": 0.992918
   }
  },
  "strength": {
   "n": 727,
   "tool": 0.284732,
   "tool_capped": 0.671252,
   "rule": 0.969739
  }
 },
 "Phenylketonuria VCEP": {
  "classifications": 323,
  "PP3": {
   "tool": {
    "evaluated": 317,
    "tp": 224,
    "fn": 30,
    "fp": 10,
    "tn": 53,
    "positive_agreement": 0.88189,
    "negative_agreement": 0.84127
   },
   "panel_rules": {
    "evaluated": 317,
    "tp": 252,
    "fn": 2,
    "fp": 41,
    "tn": 22,
    "positive_agreement": 0.992126,
    "negative_agreement": 0.349206
   }
  },
  "BP4": {
   "tool": {
    "evaluated": 10,
    "tp": 5,
    "fn": 0,
    "fp": 1,
    "tn": 4,
    "positive_agreement": 1.0,
    "negative_agreement": 0.8
   },
   "panel_rules": {
    "evaluated": 10,
    "tp": 0,
    "fn": 5,
    "fp": 0,
    "tn": 5,
    "positive_agreement": 0.0,
    "negative_agreement": 1.0
   }
  },
  "PM2": {
   "tool": {
    "evaluated": 320,
    "tp": 117,
    "fn": 179,
    "fp": 0,
    "tn": 24,
    "positive_agreement": 0.39527,
    "negative_agreement": 1.0
   },
   "panel_rules": {
    "evaluated": 320,
    "tp": 280,
    "fn": 16,
    "fp": 5,
    "tn": 19,
    "positive_agreement": 0.945946,
    "negative_agreement": 0.791667
   }
  },
  "BA1": {
   "tool": {
    "evaluated": 5,
    "tp": 0,
    "fn": 0,
    "fp": 0,
    "tn": 5,
    "positive_agreement": null,
    "negative_agreement": 1.0
   },
   "panel_rules": {
    "evaluated": 5,
    "tp": 0,
    "fn": 0,
    "fp": 0,
    "tn": 5,
    "positive_agreement": null,
    "negative_agreement": 1.0
   }
  },
  "strength": {
   "n": 321,
   "tool": 0.308411,
   "tool_capped": 0.308411,
   "rule": 0.214953
  }
 },
 "RASopathy VCEP": {
  "classifications": 233,
  "PP3": {
   "tool": {
    "evaluated": 149,
    "tp": 113,
    "fn": 0,
    "fp": 15,
    "tn": 21,
    "positive_agreement": 1.0,
    "negative_agreement": 0.583333
   },
   "panel_rules": {
    "evaluated": 149,
    "tp": 108,
    "fn": 5,
    "fp": 1,
    "tn": 35,
    "positive_agreement": 0.955752,
    "negative_agreement": 0.972222
   }
  },
  "BP4": {
   "tool": {
    "evaluated": 90,
    "tp": 21,
    "fn": 4,
    "fp": 6,
    "tn": 59,
    "positive_agreement": 0.84,
    "negative_agreement": 0.907692
   },
   "panel_rules": {
    "evaluated": 90,
    "tp": 21,
    "fn": 4,
    "fp": 1,
    "tn": 64,
    "positive_agreement": 0.84,
    "negative_agreement": 0.984615
   }
  },
  "PM2": {
   "tool": {
    "evaluated": 173,
    "tp": 90,
    "fn": 42,
    "fp": 0,
    "tn": 41,
    "positive_agreement": 0.681818,
    "negative_agreement": 1.0
   },
   "panel_rules": {
    "evaluated": 173,
    "tp": 91,
    "fn": 41,
    "fp": 0,
    "tn": 41,
    "positive_agreement": 0.689394,
    "negative_agreement": 1.0
   }
  },
  "BA1": {
   "tool": {
    "evaluated": 106,
    "tp": 0,
    "fn": 28,
    "fp": 0,
    "tn": 78,
    "positive_agreement": 0.0,
    "negative_agreement": 1.0
   },
   "panel_rules": {
    "evaluated": 106,
    "tp": 26,
    "fn": 2,
    "fp": 6,
    "tn": 72,
    "positive_agreement": 0.928571,
    "negative_agreement": 0.923077
   }
  },
  "strength": {
   "n": 167,
   "tool": 0.11976,
   "tool_capped": 0.850299,
   "rule": 0.934132
  }
 },
 "Rett and Angelman-like Disorders VCEP": {
  "classifications": 314,
  "PP3": {
   "tool": {
    "evaluated": 162,
    "tp": 72,
    "fn": 2,
    "fp": 22,
    "tn": 66,
    "positive_agreement": 0.972973,
    "negative_agreement": 0.75
   },
   "panel_rules": {
    "evaluated": 162,
    "tp": 74,
    "fn": 0,
    "fp": 5,
    "tn": 83,
    "positive_agreement": 1.0,
    "negative_agreement": 0.943182
   }
  },
  "BP4": {
   "tool": {
    "evaluated": 179,
    "tp": 57,
    "fn": 11,
    "fp": 61,
    "tn": 50,
    "positive_agreement": 0.838235,
    "negative_agreement": 0.45045
   },
   "panel_rules": {
    "evaluated": 179,
    "tp": 64,
    "fn": 4,
    "fp": 14,
    "tn": 97,
    "positive_agreement": 0.941176,
    "negative_agreement": 0.873874
   }
  },
  "PM2": {
   "tool": {
    "evaluated": 177,
    "tp": 89,
    "fn": 24,
    "fp": 1,
    "tn": 63,
    "positive_agreement": 0.787611,
    "negative_agreement": 0.984375
   },
   "panel_rules": {
    "evaluated": 177,
    "tp": 89,
    "fn": 24,
    "fp": 1,
    "tn": 63,
    "positive_agreement": 0.787611,
    "negative_agreement": 0.984375
   }
  },
  "BA1": {
   "tool": {
    "evaluated": 73,
    "tp": 0,
    "fn": 36,
    "fp": 0,
    "tn": 37,
    "positive_agreement": 0.0,
    "negative_agreement": 1.0
   },
   "panel_rules": {
    "evaluated": 73,
    "tp": 28,
    "fn": 8,
    "fp": 10,
    "tn": 27,
    "positive_agreement": 0.777778,
    "negative_agreement": 0.72973
   }
  },
  "strength": {
   "n": 260,
   "tool": 0.134615,
   "tool_capped": 0.596154,
   "rule": 0.903846
  }
 },
 "Severe Combined Immunodeficiency Disease VCEP": {
  "classifications": 313,
  "PP3": {
   "tool": {
    "evaluated": 41,
    "tp": 11,
    "fn": 0,
    "fp": 19,
    "tn": 11,
    "positive_agreement": 1.0,
    "negative_agreement": 0.366667
   },
   "panel_rules": {
    "evaluated": 41,
    "tp": 11,
    "fn": 0,
    "fp": 0,
    "tn": 30,
    "positive_agreement": 1.0,
    "negative_agreement": 1.0
   }
  },
  "BP4": {
   "tool": {
    "evaluated": 22,
    "tp": 7,
    "fn": 0,
    "fp": 3,
    "tn": 12,
    "positive_agreement": 1.0,
    "negative_agreement": 0.8
   },
   "panel_rules": {
    "evaluated": 22,
    "tp": 7,
    "fn": 0,
    "fp": 0,
    "tn": 15,
    "positive_agreement": 1.0,
    "negative_agreement": 1.0
   }
  },
  "PM2": {
   "tool": {
    "evaluated": 245,
    "tp": 48,
    "fn": 146,
    "fp": 0,
    "tn": 51,
    "positive_agreement": 0.247423,
    "negative_agreement": 1.0
   },
   "panel_rules": {
    "evaluated": 245,
    "tp": 141,
    "fn": 53,
    "fp": 2,
    "tn": 49,
    "positive_agreement": 0.726804,
    "negative_agreement": 0.960784
   }
  },
  "BA1": {
   "tool": {
    "evaluated": 79,
    "tp": 15,
    "fn": 22,
    "fp": 0,
    "tn": 42,
    "positive_agreement": 0.405405,
    "negative_agreement": 1.0
   },
   "panel_rules": {
    "evaluated": 79,
    "tp": 37,
    "fn": 0,
    "fp": 0,
    "tn": 42,
    "positive_agreement": 1.0,
    "negative_agreement": 1.0
   }
  },
  "strength": {
   "n": 50,
   "tool": 0.14,
   "tool_capped": 0.9,
   "rule": 1.0
  }
 }
}
```
