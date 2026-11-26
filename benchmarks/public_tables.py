"""Write the per-variant result tables that can be redistributed (benchmarks/results/public/).

    python benchmarks/public_tables.py

The benchmark scripts write their per-variant tables to benchmarks/results/ with every score they used,
including scores whose terms do not allow redistribution under this repository's GPL-3.0: REVEL, CADD and
SpliceAI (non-commercial use only), the dbNSFP predictors (several for academic use only), EVE and popEVE,
COSMIC and GENIE counts, and the ProteinGym assay values (taken from the original studies). Those full
tables stay local (git-ignored); anyone can regenerate them by running the scripts with their own access
to the sources.

This script keeps, for each table, the variant keys, labels and dates, MissensExplorer's own outputs
(scores, calls, criteria, points), and the scores that may be redistributed with attribution:
AlphaMissense (CC BY 4.0; Cheng et al., Science 2023), gnomAD frequencies (gnomAD's data are free to
reuse) and ClinGen/ClinVar data (CC0 / public domain). Every other score column is dropped, and a NOTICE
lists the sources and their terms.
"""

from __future__ import annotations

import re
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))

from compare_tools import TOOLS

RESULTS = Path(__file__).resolve().parent / "results"
PUBLIC = RESULTS / "public"
KEEP_SCORES = {"AlphaMissense", "alphamissense"}
THIRD_PARTY = (
    {name for name in TOOLS if name not in KEEP_SCORES}
    | {"revel", "cadd_phred", "spliceai_max", "eve", "popeve", "sift_score", "polyphen_score"}
    | {"cosmic", "genie", "dms", "bin"}
)
# third-party scores quoted inside text columns ("... SpliceAI 0.15 does not rule out ...")
SCORE_IN_TEXT = re.compile(r"\b(SpliceAI|REVEL|CADD|EVE|popEVE|BayesDel|MetaRNN|ClinPred)( score)?\s+-?\d+(?:\.\d+)?")
TABLES = (
    "clinvar_benchmark_variants.csv.gz",
    "clinvar_oncogenicity_variants.csv.gz",
    "egfr_validation.csv",
    "egfr_oncogenicity_variants.csv.gz",
    "erepo_validation_variants.csv.gz",
    "erepo_vcep_rules_variants.csv.gz",
    "erepo_clinvar_history_variants.csv.gz",
    "model_training_set.csv.gz",
    "proteingym_variants.csv.gz",
    "tool_comparison_variants.csv.gz",
)
NOTICE = """Per-variant results of the MissensExplorer benchmarks (redistributable subset)

Written by benchmarks/public_tables.py from the full tables that the benchmark scripts produce. Columns
with third-party scores whose terms do not allow redistribution under GPL-3.0 were dropped: REVEL, CADD,
SpliceAI, EVE, popEVE, SIFT, PolyPhen-2 and every dbNSFP predictor, COSMIC/GENIE counts and ProteinGym
assay values and bins; scores quoted in text columns are replaced by "[score removed]". Regenerate the
full tables by running the scripts (see README.md).

Kept, with their sources and terms:
- variant keys, labels and dates from ClinVar (public domain) and the ClinGen Evidence Repository (CC0);
- MissensExplorer's own scores, calls, ACMG/AMP and oncogenicity criteria and points (GPL-3.0 project);
- AlphaMissense scores (CC BY 4.0; Cheng et al., Science 2023);
- gnomAD allele frequencies (gnomAD, Broad Institute; free to reuse);
- ClinGen CSpec-derived codes (CC0).
"""


def main() -> int:
    PUBLIC.mkdir(parents=True, exist_ok=True)
    for name in TABLES:
        src = RESULTS / name
        if not src.exists():
            print(f"  missing: {name}")
            continue
        df = pd.read_csv(src, dtype=str, keep_default_na=False)
        drop = [c for c in df.columns if c in THIRD_PARTY]
        out = df.drop(columns=drop)
        for c in out.columns:  # scores quoted in free text
            out[c] = out[c].str.replace(SCORE_IN_TEXT, r"\1 [score removed]", regex=True)
        dest = PUBLIC / (name if name.endswith(".gz") else name + ".gz")
        out.to_csv(dest, index=False, compression="gzip")
        print(f"  {name}: {len(out):,} rows, {len(out.columns)} columns kept, {len(drop)} dropped")
    (PUBLIC / "NOTICE.txt").write_text(NOTICE)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
