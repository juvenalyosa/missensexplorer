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
- Desktop app: the welcome page describes the prediction as the model (AlphaMissense, REVEL, CADD
  and popEVE) instead of the weighted vote it replaced in 2.2.
- Desktop and web apps: variants without PP3/BP4 (everything that is not missense, such as a
  stop-gained variant) no longer show a note about REVEL's PP3/BP4 thresholds under their ACMG
  evidence.
- `missensexplorer cache --help` prints the help. Before, it printed the cache path, and any other
  word was taken as "path". `missensexplorer cache path` is now also accepted.
- `missensexplorer web` sets the upload limit (1 GB), the toolbar and the theme colours itself.
  Before, they applied only when the app was started from the repository folder, and elsewhere
  Streamlit's 200 MB upload limit applied.

- A protein substitution that several nucleotide changes produce, such as KRAS Q61H (CAA to CAC or
  CAT), used to be rejected by Ensembl. It is now asked again with `ambiguous_hgvs` and gives one row
  per nucleotide change, with a note.
- Insertions and deletions typed as HGVS (`CFTR:c.1521_1523del`) are re-annotated in their
  left-aligned VCF form. Before, the position Ensembl resolved them to missed CADD's indel scores,
  and sometimes gnomAD's frequency and the ClinVar record, that the same variant gets from a VCF.
  An insertion or deletion written as a protein change (`p.Phe508del`) now says that other DNA
  changes can give the same protein and which one Ensembl mapped.
- GRCh37: UniProt features and the AlphaFold structure were missing, because Ensembl's GRCh37 server
  reports UniProt entry names (BRAF_HUMAN); they are now converted to accessions. PS1/PM5 said
  "none" for rows without a MANE Select transcript (all of GRCh37), where nothing could be compared;
  it now says "not checked" and why.
- VCF records removed by *FILTER = PASS only* are now counted in the skipped note of the desktop and
  web apps and in the terminal summary, which also lists the other skipped records.
- Ensembl resolves a protein change on any isoform of the gene, so a change written in another
  numbering (MUTYH p.Gly396Asp) is reported at its position on the MANE transcript (G368D). Such rows
  now say so, and say *check the input* when the amino acids differ too (TP53 p.Ala175His is C182S on
  the MANE transcript), instead of changing the variant silently.
- Readable messages replace Ensembl's internal errors for a position outside the protein, a wrong
  reference amino acid and a frameshift written as a protein change. Errors cached by older versions
  are shown with the new wording.

- A protein change given with a gene symbol is resolved on the gene's canonical (MANE Select)
  transcript. Ensembl resolved it on the first isoform with a matching residue, which can number the
  protein differently. `EGFR T725M`, `EGFR L747P` and `EGFR p.Glu746_Ala750del` landed in the 3' UTR and
  gave no row and no error; `EGFR V689L` became V742L. The change is now asked again on the canonical
  transcript when VEP's first answer does not carry it there; the first answer stays when the canonical
  transcript has another residue at that position (MUTYH p.Gly396Asp, with its note).
- A protein change that still gives no protein-coding row is reported as an error instead of being
  dropped silently.
- Protein insertions and deletion-insertions (`EGFR p.Asp770_Asn771insSerValAsp`), which Ensembl cannot
  trace back to DNA, say so and ask for the c. change. Before, they said the reference amino acid was
  wrong, or that VEP could not map them.
- A protein deletion that left-alignment moves into the previous codon (`EGFR p.Leu747del`) no longer
  gets the note about another isoform's numbering.
- VEP's own error message is kept in the cache, so a later version can reword it.

### Added
- EGFR validation (`benchmarks/validate_egfr.py`): 356 tumour variants against a laboratory
  transformation assay (Ikeuchi et al., Oncogene 2022). The scores detect loss of function, not
  activation; see `docs/model.md`. OncoKB's annotation from the paper's table is used only for a
  summary AUC and is not written out per variant, because OncoKB's terms do not allow sharing it.
- Somatic oncogenicity points, following the ClinGen/CGC/VICC standard (Horak et al., Genet Med
  2022). They are optional: `--somatic` on the command line, or a switch in the desktop and web apps.
  - Only the criteria that can be computed automatically are applied, with open data:
    - hotspots from cancerhotspots.org;
    - gene roles from IntOGen (CC0);
    - functional studies curated in CIViC (CC0) for OS2/SBS2;
    - ClinVar oncogenicity classifications of other variants at the residue;
    - gnomAD and the model.
  - The model's call counts (OP1/SBP1) only in genes that are not activating. In EGFR it pointed the
    wrong way.
  - CIViC's clinical evidence is shown next to the points but adds none: the best evidence level (A
    to E), drugs and tumour type.
  - New columns: `oncogenicity`, `oncogenicity_points`, `oncogenicity_codes`,
    `oncogenicity_evidence`, `clinical_evidence_level` and `clinical_evidence`.
  - Validated on EGFR and TP53 laboratory data (`benchmarks/validate_oncogenicity.py`).
- ClinVar records keep the review status of their oncogenicity classification. The cache key changed,
  so records are fetched again once.
- A tutorial [wiki](https://github.com/juvenalyosa/missensexplorer/wiki) with screenshots of every
  feature, and worked examples of every input format and kind of variant.

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
