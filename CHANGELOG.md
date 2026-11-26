# Changelog

The format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/). Versions follow
[Semantic Versioning](https://semver.org/).

## [2.5.0] — 2026-11-26

Corrections prompted by six rounds of review of the evaluation. The model and its PP3/BP4
thresholds are retrained and recalibrated, so the results of 2.4 change.

### Fixed
- Missense model: the genes of the ClinVar benchmark were left out of training by an exact match on
  ClinVar's gene-symbol field, which lists overlapping genes together ('LDLR;MIR6886'). Training now
  leaves a variant out when the gene of its ClinVar name or any of its symbols is a benchmark gene;
  cross-validation and the calibration's gene folds use the gene of the ClinVar name. The model is
  retrained (version 2.5) and the PP3/BP4 thresholds recalibrated.
- A missing predictor score is filled with the training mean, so it adds nothing to the score. The
  median used before added a small term that the explanation did not show.
- BP4 is withheld when SpliceAI does not rule out a splicing effect (score above 0.1, the ClinGen
  splicing subgroup's threshold for no impact): BP4 means no predicted impact on the protein and on
  splicing.
- Evidence worth 3 points is written `PP3_Moderate+` / `BP4_Moderate+` instead of `_Moderate`, which
  a curator would read as 2 points.
- BA1 needs a population of at least 2,000 alleles (ClinGen, Ghosh et al. 2018): gnomAD itself is
  asked about variants at or above 5%, and BA1 is not given when no such population exceeds 5%. BA1
  now needs a frequency above 5%, as the guideline says (it was at or above).
- BP4 is also withheld when there is no SpliceAI score to rule out a splicing effect; PS1/PM5 say so.
- BA1 is not given below the higher BA1 threshold of a gene's expert panels (VWF 10%, ABCA4 16.3%):
  BA1 at 5% applies to genes without a gene-specific rule (Ghosh et al. 2018); nor does the prediction
  column call such a variant "common". When gnomAD could not be asked for allele numbers, the BA1 text
  says so.
- A truncation that spans two exons ("9-10/10" in VEP) crosses an exon-intron boundary: its
  nonsense-mediated decay is "uncertain" and OVS1 is graded by the share of the protein removed (before,
  its decay was not assessed and it got full OVS1).
- When Ensembl cannot be reached but the cache holds earlier answers (a rerun of a published analysis),
  the run uses every cached answer, whatever its age, and says so (`cache_only` in the run record),
  instead of stopping; this also holds when Ensembl answers with server errors or times out.
- A frequent variant that ClinVar classifies as pathogenic with review (two stars or more), such as HbS
  in HBB, carries a caution in `acmg_notes`: panels keep such variants pathogenic on other evidence, so
  BA1 and BS1 need review. The classification is shown, never counted as evidence.
- PS1/PM5: a ClinVar reference variant with a single review star counts one step lower (a pathogenic
  one like a likely pathogenic one; a likely pathogenic one not at all). Expert panels agreed with PM5
  from single-star references in only about a third of their classifications.
- Somatic oncogenicity:
  - gene roles: IntOGen's role weighted by mutated samples in cohorts with a role (the majority of
    cohorts made SMARCA4, KEAP1, GATA3 and CHEK2 "activating"), when it rests on at least 2 cohorts with
    a role or 20 samples in them, and never on fewer than 10 samples (POLE was "activating" on 4, CDKN2C
    on 8), not outnumbered by samples in cohorts of ambiguous role (which made PMS2 and SDHA "activating";
    cohorts of ambiguous role no longer count towards the two cohorts, which made MAML2 loss-of-function
    on one cohort); UniProt's curated keyword (tumour suppressor or proto-oncogene) when
    IntOGen is absent, ambiguous or too thin; when they disagree, IntOGen decides only with at least 3
    cohorts and 50 samples, otherwise the role is ambiguous (UniProt's Proto-oncogene keyword also marks
    fusion genes such as RUNX1, NSD1 and CBFB, whose point mutations inactivate them; a few IntOGen
    samples made SMARCB1 "activating"); genes known only to UniProt count too, and a note says when
    UniProt could not be reached;
  - OS1 never lowers the points: a ClinVar classification of the same amino-acid change (with OM3 or
    OP3, which the standard combines with OS1) and the evidence the standard does not use with OS1 (OS3,
    OS2, OM4, with OM3 or OP3) are alternatives, and the stronger is kept (before, a likely oncogenic
    reference, OS1_Moderate, removed OS3, OS2 and OM4: TP53 p.K132N fell to VUS; and OS1 removed OM3 and
    OP3). OM3 comes back when OS1 replaces OM4. OS1_Moderate, the tool's own level for a likely oncogenic
    reference, is not combined with OM3 or OP3 (that reference may itself rest on the same hotspot) and
    gives way to the standard's criteria on a tie. Keeping the excluded evidence when it is stronger
    departs from the standard; the evidence then gives the points under the standard's rule;
  - OP4 ("absent from gnomAD") needs gnomAD itself to confirm the absence, with enough coverage, as
    PM2 does (before, it was also given when gnomAD had not been asked);
  - CIViC is read from its monthly release (citable), not the nightly file, with English month names
    whatever the locale (the desktop app in Spanish built wrong release names);
  - OVS1 is graded as ClinGen's loss-of-function guidance grades PVS1, which the standard
    recommends: truncations that escape nonsense-mediated decay get OVS1_Strong or OVS1_Moderate by
    the share of protein removed (before: no OVS1 anywhere in the last exon, which missed most APC
    truncations); canonical splice variants by whether skipping the exon keeps the frame; lost start
    codons OVS1_Moderate;
  - OP1/SBP1 are applied in every gene, as the standard does, from the calibrated PP3/BP4 evidence
    (before: not in activating genes, and from the model's call cut-offs);
  - OS2 uses only CIViC's preclinical (functional) evidence, not clinical observations;
  - cancerhotspots.org gene symbols are updated to HGNC's current ones (H3F3A is H3-3A), so those
    hotspots are found;
  - a variant's own ClinVar oncogenicity classification is found by its position, not only through
    VEP's identifiers, so it is shown for many more variants (it is still never counted as evidence).
- Desktop app: the *With frequency* and *AlphaMissense* column headers are no longer cut off.

### Added
- ClinGen expert panel specifications (`vcep.py`; table built from ClinGen's Criteria Specification
  Registry by `benchmarks/make_vcep_table.py`: 127 released specifications, 126 with genes, 138
  genes). For a variant in a covered gene, the new `vcep_specification` column names every
  specification with its disease and link, and the new `acmg_notes` column gives one note per
  specification on what differs from the tool's evidence:
  - a strength the panel does not allow (stronger, or the 3-point level, which no panel uses; a weaker
    PM2; a PS1/PM5 the tool lowered for a likely pathogenic or single-star reference, for which the
    panel's specification has no level);
  - PP3/BP4 that the panel does not use for missense variants although CSpec marks them applicable
    (the SCID genes, PALB2, APC, VHL, DDC, CDH1 and the brain malformation genes), quoted;
  - that the panel sets its own predictor and thresholds for PP3/BP4 (not for panels that use neither);
  - the panel's own BA1 threshold, read from its text (VWF's is above 10%, most are far below 5%), with a
    warning when the panel applies it to the filtering allele frequency, which the tool does not compute,
    and the text quoted when it sets further conditions (allele counts, males only, other datasets,
    excluded or founder variants), whichever frequency the panel uses; a variant on the panel's own BA1
    exclusion list (RYR1, NEB, ACTA1; matched on MANE Select, NEB's list mapped from NM_001271208.2, and
    shown with its transcript) or on ClinGen's BA1 exception list is said not to meet BA1 (nor BS1 where
    the panel says so);
  - its conditions the tool cannot check, quoted from the clause that states them (for example,
    references classified under its specification, PM5 not with PM1);
  - the disease and mode of inheritance of each specification, from CSpec's structured fields;
  - for the ten panels whose PP3/BP4, PM2 and BA1 rules were encoded by hand (`vcep_rules.py`), what
    those rules give for the variant (PM2 "not assessed" without gnomAD data; the Myeloid Malignancy
    panel's 20x coverage for an absent variant).
  The notes are shown in the ACMG section of the desktop app, the web app and the reports.
- A note when PP3 at moderate or stronger meets PS1/PM5, which may share residue-level information
  (Pejaver et al. 2022); it is in `acmg_notes` too.
- The run record lists the releases of the data sets used (CIViC, UniProt, IntOGen, cancerhotspots.org,
  CSpec); the somatic benchmark results keep the full record of every run.
- Benchmarks: per-variant tables for the ClinVar benchmark, the tool comparison and ProteinGym; the
  ProteinGym comparison is paired by variant and reported by assay type; the tool comparison adds
  variants absent from gnomAD and below 0.01%; the calibration tests leave out benchmark variants that
  are in the calibration set and give exact (Katz) intervals and the probability that each band
  reaches its ratio; the eRepo comparison (`validate_erepo.py`, `erepo_analyses.py`) and the ClinVar
  oncogenicity comparison (`validate_clinvar_oncogenicity.py`) are new; `erepo_vcep_rules.py` applies the
  hand-encoded rules of the ten largest panels to their classifications (with the tool capped at each
  specification's maximum strength, by approval period, with intervals over panels), and
  `erepo_clinvar_history.py` recomputes PS1/PM5 with the monthly ClinVar release before each panel's
  decision plus the panels' own earlier classifications. `public_tables.py` writes the per-variant
  tables without the third-party scores that may not be redistributed (`benchmarks/results/public/`);
  the full tables are git-ignored (the two that 2.4 tracked, `model_training_set.csv.gz` and
  `egfr_validation.csv`, leave the index with this release); scores and cancerhotspots.org tumour counts
  quoted in text columns are removed, `NOTICE.txt` lists each kept source and its terms, and the files are
  byte-for-byte reproducible. `benchmarks/requirements-lock.txt` lists the exact package versions of the
  analysis environment, and the `benchmark` extra now includes SciPy, which the scripts need.
- CI also runs the desktop app's tests and the CIViC test under a Spanish locale.
- The hotspot-only comparison of the EGFR/TP53 evaluation takes the hotspot criteria met from the run without
  ClinVar (a reference could remove them from the final points); the ClinVar oncogenicity evaluation gives
  intervals from a bootstrap of genes and keeps genes whose IntOGen role rests on too few samples apart.

### Changed
- The `acmg_notes` column is new and sits after `vcep_specification`; the overlap and expert-panel
  notes moved there from `notes`. Scripts that read the CSV by column position must be updated.

## [2.4.0] — 2026-11-25

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
- Double-click installers with everything inside (Python, Qt, the 3D viewer). Run
  `packaging/build.py` on each system, or the `installers` workflow, which builds all of them and
  self-tests each package:
  - macOS: a .dmg for Apple Silicon and one for Intel (macOS 13 or later);
  - Windows: a Setup.exe that installs for the user, with no admin rights;
  - Linux: a .deb for Ubuntu/Debian and a .tar.gz with an install script.

  The installers pin Qt 6.11. Qt 6.8 to 6.10 would support older systems, but their web engine
  crashes on current macOS.
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
