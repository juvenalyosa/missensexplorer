# PS1/PM5 against the expert panels, with the references of each panel's decision

Produced by `benchmarks/erepo_clinvar_history.py`. `tool`: the tool's calls (current ClinVar, MANE transcript); `now`: this re-implementation on the 22 Nov 2026 release (panel's transcript); `then`: the same rule with the last monthly release before the panel's approval month plus the panels' own classifications approved by that day; `then_clinvar`: the archived release alone.

```
{
 "variants": 6310,
 "transcripts": 145,
 "releases_archived": 106,
 "releases_used": {
  "2016-11": 23,
  "2017-03": 60,
  "2017-04": 9,
  "2017-09": 4,
  "2018-03": 1,
  "2018-06": 4,
  "2018-07": 39,
  "2018-08": 15,
  "2018-09": 4,
  "2018-10": 2,
  "2018-11": 39,
  "2018-12": 2,
  "2019-01": 4,
  "2019-02": 4,
  "2019-03": 15,
  "2019-04": 18,
  "2019-05": 10,
  "2019-06": 13,
  "2019-07": 7,
  "2019-08": 6,
  "2019-09": 2,
  "2019-10": 9,
  "2019-11": 14,
  "2019-12": 12,
  "2020-01": 10,
  "2020-02": 20,
  "2020-03": 9,
  "2020-04": 15,
  "2020-05": 31,
  "2020-06": 22,
  "2020-07": 15,
  "2020-08": 8,
  "2020-09": 10,
  "2020-10": 10,
  "2020-11": 2,
  "2020-12": 10,
  "2021-01": 11,
  "2021-02": 76,
  "2021-03": 8,
  "2021-04": 17,
  "2021-05": 54,
  "2021-06": 9,
  "2021-07": 10,
  "2021-08": 20,
  "2021-09": 20,
  "2021-10": 12,
  "2021-11": 39,
  "2021-12": 25,
  "2022-01": 40,
  "2022-02": 46,
  "2022-03": 31,
  "2022-04": 33,
  "2022-05": 23,
  "2022-06": 25,
  "2022-07": 44,
  "2022-08": 17,
  "2022-09": 33,
  "2022-10": 4,
  "2022-11": 34,
  "2022-12": 16,
  "2023-01": 44,
  "2023-02": 30,
  "2023-03": 250,
  "2023-04": 22,
  "2023-05": 63,
  "2023-06": 25,
  "2023-07": 145,
  "2023-08": 32,
  "2023-09": 40,
  "2023-10": 46,
  "2023-11": 68,
  "2023-12": 121,
  "2024-01": 142,
  "2024-02": 61,
  "2024-03": 114,
  "2024-04": 93,
  "2024-05": 97,
  "2024-06": 146,
  "2024-07": 253,
  "2024-08": 264,
  "2024-09": 156,
  "2024-10": 86,
  "2024-11": 120,
  "2024-12": 279,
  "2025-01": 118,
  "2025-02": 83,
  "2025-03": 74,
  "2025-04": 188,
  "2025-05": 112,
  "2025-06": 101,
  "2025-07": 87,
  "2025-08": 125,
  "2025-09": 118,
  "2025-10": 174,
  "2025-11": 173,
  "2025-12": 132,
  "2026-01": 84,
  "2026-02": 56,
  "2026-03": 149,
  "2026-04": 109,
  "2026-05": 108,
  "2026-06": 73,
  "2026-07": 115,
  "2026-08": 137,
  "2026-09": 93,
  "2026-10": 79
 },
 "reimplementation_agrees_with_tool": 6283,
 "reimplementation_agrees_with_tool_where_applied": [
  1774,
  1801
 ],
 "tool_calls": {
  "PM5_Supporting": 676,
  "PM5": 982,
  "PS1": 74,
  "PS1_Moderate": 61
 },
 "calls_now": {
  "PM5_Supporting": 675,
  "PM5": 975,
  "PS1": 72,
  "PS1_Moderate": 61
 },
 "calls_at_approval": {
  "PM5_Supporting": 510,
  "PM5": 881,
  "PS1": 57,
  "PS1_Moderate": 49
 },
 "tool_call_without_reference_at_approval": 354,
 "tool_calls_total": 1793,
 "agreement": {
  "tool": {
   "PS1": {
    "tp": 39,
    "fn": 26,
    "fp": 35,
    "tn": 1960,
    "positive_agreement": 0.6,
    "negative_agreement": 0.982456,
    "ppv": 0.527027
   },
   "PM5": {
    "tp": 641,
    "fn": 83,
    "fp": 527,
    "tn": 1994,
    "positive_agreement": 0.885359,
    "negative_agreement": 0.790956,
    "ppv": 0.548801
   }
  },
  "now": {
   "PS1": {
    "tp": 39,
    "fn": 26,
    "fp": 34,
    "tn": 1961,
    "positive_agreement": 0.6,
    "negative_agreement": 0.982957,
    "ppv": 0.534247
   },
   "PM5": {
    "tp": 637,
    "fn": 87,
    "fp": 523,
    "tn": 1998,
    "positive_agreement": 0.879834,
    "negative_agreement": 0.792543,
    "ppv": 0.549138
   }
  },
  "then": {
   "PS1": {
    "tp": 35,
    "fn": 30,
    "fp": 32,
    "tn": 1963,
    "positive_agreement": 0.538462,
    "negative_agreement": 0.98396,
    "ppv": 0.522388
   },
   "PM5": {
    "tp": 606,
    "fn": 118,
    "fp": 413,
    "tn": 2108,
    "positive_agreement": 0.837017,
    "negative_agreement": 0.836176,
    "ppv": 0.594701
   }
  },
  "then_strict": {
   "PS1": {
    "tp": 29,
    "fn": 36,
    "fp": 28,
    "tn": 1967,
    "positive_agreement": 0.446154,
    "negative_agreement": 0.985965,
    "ppv": 0.508772
   },
   "PM5": {
    "tp": 538,
    "fn": 186,
    "fp": 350,
    "tn": 2171,
    "positive_agreement": 0.743094,
    "negative_agreement": 0.861166,
    "ppv": 0.605856
   }
  },
  "then_clinvar": {
   "PS1": {
    "tp": 29,
    "fn": 36,
    "fp": 27,
    "tn": 1968,
    "positive_agreement": 0.446154,
    "negative_agreement": 0.986466,
    "ppv": 0.517857
   },
   "PM5": {
    "tp": 508,
    "fn": 216,
    "fp": 337,
    "tn": 2184,
    "positive_agreement": 0.701657,
    "negative_agreement": 0.866323,
    "ppv": 0.601183
   }
  }
 },
 "pm5_paired": {
  "panel_pm5_met": 724,
  "lost_then": 51,
  "gained_then": 20,
  "mcnemar_p_unclustered": 0.000303163413558493,
  "positive_agreement_change_ci": [
   -0.0803,
   0.0028
  ],
  "ppv_change_ci": [
   0.0134,
   0.0887
  ]
 },
 "pm5_attribution": {
  "tool_pm5_against_panel": 527,
  "gross_no_pm5_reference_then": 138,
  "gross_share": 0.2619,
  "discordant_now": 523,
  "discordant_removed_then": 132,
  "concordant_now": 637,
  "concordant_removed_then": 51,
  "excess_share": 0.1723,
  "clinvar_only_discordant_removed": 207,
  "clinvar_only_concordant_removed": 143,
  "reimplementation_basis": {
   "discordant_lost": [
    132,
    523
   ],
   "concordant_lost": [
    51,
    637
   ],
   "gross": 0.2524,
   "net": 0.1723,
   "gross_ci": [
    0.1578,
    0.3674
   ],
   "net_ci": [
    0.0684,
    0.3002
   ],
   "panels": 47,
   "strictly_earlier_panel_references": {
    "discordant_lost": [
     195,
     523
    ],
    "concordant_lost": [
     118,
     637
    ],
    "gross": 0.3728,
    "net": 0.1876
   }
  },
  "by_period": {
   "2016-2021": {
    "discordant_lost": [
     60,
     127
    ],
    "concordant_lost": [
     2,
     79
    ],
    "gross": 0.4724,
    "net": 0.4471,
    "gross_ci": [
     0.2333,
     0.6023
    ],
    "net_ci": [
     0.1799,
     0.5741
    ],
    "panels": 14
   },
   "2022-2023": {
    "discordant_lost": [
     30,
     135
    ],
    "concordant_lost": [
     14,
     127
    ],
    "gross": 0.2222,
    "net": 0.112,
    "gross_ci": [
     0.1346,
     0.5
    ],
    "net_ci": [
     0.0279,
     0.3975
    ],
    "panels": 23
   },
   "2024-2026": {
    "discordant_lost": [
     42,
     261
    ],
    "concordant_lost": [
     35,
     431
    ],
    "gross": 0.1609,
    "net": 0.0797,
    "gross_ci": [
     0.1103,
     0.2262
    ],
    "net_ci": [
     0.0189,
     0.1544
    ],
    "panels": 45
   },
   "2016-2021 without Phenylketonuria": {
    "discordant_lost": [
     20,
     58
    ],
    "concordant_lost": [
     2,
     42
    ],
    "gross": 0.3448,
    "net": 0.2972,
    "gross_ci": [
     0.2021,
     0.65
    ],
    "net_ci": [
     0.1429,
     0.5806
    ],
    "panels": 13
   }
  },
  "phenylketonuria_2016_2021": {
   "discordant_lost": [
    40,
    69
   ],
   "concordant_lost": [
    0,
    37
   ],
   "gross": 0.5797,
   "net": 0.5797
  },
  "within_panel_net": {
   "net": 0.2067,
   "discordant_in_panels_used": [
    485,
    523
   ]
  },
  "release_lag_days": {
   "median": 46.0,
   "min": 28,
   "max": 61
  },
  "clinvar_only_excess_share": 0.1713
 },
 "without_transcript_switches": {
  "flagged_classifications": 47,
  "flagged_without_reference_then": 10,
  "then": {
   "PS1": {
    "tp": 35,
    "fn": 30,
    "fp": 32,
    "tn": 1945,
    "positive_agreement": 0.538462,
    "negative_agreement": 0.983814,
    "ppv": 0.522388
   },
   "PM5": {
    "tp": 606,
    "fn": 112,
    "fp": 411,
    "tn": 2086,
    "positive_agreement": 0.844011,
    "negative_agreement": 0.835402,
    "ppv": 0.59587
   }
  }
 }
}
```
