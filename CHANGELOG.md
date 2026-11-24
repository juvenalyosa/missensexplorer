# Changelog

The format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/). Versions follow
[Semantic Versioning](https://semver.org/).

## [Unreleased]

### Changed
- Code-quality pass with no change in results:
  - every module, class and function is typed and documented (Google-style docstrings);
  - the long functions of the pipeline and the desktop app are split into named steps;
  - `ruff` (lint and formatting) and `mypy` run in CI and as pre-commit hooks;
  - `CONTRIBUTING.md`, `docs/architecture.md` and this changelog are new.

### Fixed
- A service that keeps failing, such as AlphaFold DB rate-limiting with HTTP 429, no longer stalls a
  run. After three consecutive failed requests the host is skipped for the rest of the run, and each
  request waits at most 90 s in total. Before, a 200-protein run could wait more than half an hour.
  There is also no wait after a request's last attempt.
- Each analysis closes its cache connection. Before, the connections stayed open until garbage
  collection.
- `classify` treats NaN scores as missing, so results-table rows no longer give a basis such as
  "CADD nan".
- A warning is no longer repeated when the genome build is switched automatically with a gene
  filter.
- Desktop app: the *Prediction* column sorts consistently when a variant has no score. The popEVE
  caption now goes through the translation table.
- A VEP error response that is JSON but not an object is now read as an error message.

## [2.3.0] — 2026-11-23

### Changed
- The missense model is retrained on a set enriched in rare benign variants (9,218 variants, of
  which 1,826 are rare and benign), because v2.2 called more rare benign variants damaging than the
  v2.1 vote.
- The damaging cut-off is now 0.76, chosen by cross-validation split by gene. Variants between 0.70
  and 0.76 are now *uncertain*, never *benign*.

### Added
- Comparison with 31 other missense predictors on the held-out ClinVar benchmark
  (`benchmarks/compare_tools.py`). The model is the best of the tools that do not use allele
  frequency: AUC 0.985.
- Independent validation on deep mutational scanning data (ProteinGym, `benchmarks/validate_proteingym.py`):
  first of 34 scores. Inputs tested and rejected are documented in `benchmarks/feature_tests.py`.
- ACMG PP3/BP4 from the model score, with thresholds calibrated by ClinGen's method
  (`benchmarks/calibrate_acmg.py`). `--acmg-source revel` keeps ClinGen's published REVEL thresholds.
- ACMG PS1/PM5 from ClinVar's pathogenic variants at the same residue.
- PM2_Supporting from a direct gnomAD check of absence and coverage (`--no-gnomad-check` skips it).
- `acmg_codes` summarises the criteria met, with their points.
- `score_with_frequency`: the model score combined with gnomAD grpmax, for ranking variants (AUC
  0.996). It is not ACMG evidence.

## [2.2.0] — 2026-11-22

### Added
- Missense verdict from a logistic regression over AlphaMissense, REVEL, CADD and popEVE, trained on
  expert-reviewed ClinVar variants with every benchmark gene held out. Every call shows each
  predictor's contribution.
- ACMG PP3/BP4 with ClinGen-calibrated intervals for REVEL (Pejaver et al. 2022), else AlphaMissense
  (Bergquist et al. 2025).
- Population-frequency rule (BA1 at 5% gnomAD grpmax), with ClinGen's BA1 exception list.

## [2.1.0] — 2026-11-22

### Added
- PySide6 desktop app (English and Spanish) and macOS app packaging.
- Held-out ClinVar benchmark (`benchmarks/clinvar_benchmark.py`).
- Per-source status, run record (versions, options, input checksum), local VEP backend and
  streaming VCF reading.
- Recorded edge-case fixtures, and nightly contract tests against the live APIs.

### Changed
- SIFT and PolyPhen-2 vote only when no modern predictor has a score; the label then says *limited
  evidence*.

## [2.0.0] — 2026-11-22

### Added
- Rewrite around batched Ensembl VEP REST calls. Modern predictors (AlphaMissense, REVEL, CADD,
  EVE/popEVE, SpliceAI), exact-allele ClinVar, UniProt features, AlphaFold pLDDT and mutant protein
  sequences.
- Inputs: free text (protein changes, HGVS, rsIDs, genomic positions), VCF/VCF.gz and CSV/TSV/Excel,
  with automatic genome-build detection.
- Command line and Streamlit web app.

### Fixed
- Several bugs of the 1.x scripts (kept in `legacy/`), listed in the README under *What changed in
  version 2*.

## [1.0] — 2024-04-30

- Original scripts (UCSC, Ensembl, UniProt and NCBI web requests per SNP), now in `legacy/`.
