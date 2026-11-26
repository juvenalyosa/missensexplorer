# Somatic oncogenicity points on laboratory data (EGFR and TP53)

Produced by `benchmarks/validate_oncogenicity.py`. The points follow the ClinGen/CGC/VICC standard (Horak et al., Genet Med 2022) with the criteria computed automatically (see `src/missensexplorer/oncogenicity.py`): hotspots from cancerhotspots.org, gene roles from IntOGen (IntOGen-Drivers-20240920) or UniProt keywords, ClinVar oncogenicity classifications of other variants at the residue, functional evidence curated in CIViC (OS2/SBS2), gnomAD (OP4 only when gnomAD itself confirms the absence), and the calibrated computational evidence (OP1/SBP1, in every gene).

- **EGFR** (oncogene): Ikeuchi et al., Oncogene 2022 (CC BY 4.0): NIH-3T3 transformation (3-4 vs 1) and Ba/F3 transformation.
- **TP53** (tumour suppressor): NCI TP53 Database r21: Kato et al. 2003 yeast transactivation (non-functional vs functional/supertrans) and Giacomelli et al. 2018 loss of function.
- **AUC**: chance that a variant positive in the assay has more points than a negative one (ties count half), with a 95% bootstrap interval. *Without ClinVar*: the same run with ClinVar off.
- The points are an integer scale: many variants share a value, which caps the AUC.

## Ranking

| Gene | Assay | Variants | Positives | Model alone | Points | Points without ClinVar |
|---|---|---|---|---|---|---|
| EGFR | 3T3 (missense) | 225 | 38 | 0.542 | **0.693** (0.581–0.797) | 0.697 |
| EGFR | 3T3 (all) | 243 | 46 | 0.542 | **0.723** (0.619–0.813) | 0.726 |
| EGFR | Ba/F3 (missense) | 258 | 200 | 0.214 | **0.396** (0.341–0.45) | 0.395 |
| EGFR | Ba/F3 (all) | 279 | 212 | 0.214 | **0.393** (0.335–0.454) | 0.39 |
| TP53 | Kato (missense) | 1724 | 493 | 0.927 | **0.931** (0.916–0.947) | 0.928 |
| TP53 | Giacomelli (missense) | 2158 | 570 | 0.84 | **0.848** (0.826–0.868) | 0.845 |

## CIViC

Items that cite the validation assays or yeast assays like Kato's are left out here (not in the app): PMID 35304574, Ikeuchi et al. 2022 (the EGFR assay); PMID 12826609, Kato et al. 2003 (the TP53 yeast assay); PMID 30224644, Giacomelli et al. 2018 (the TP53 human-cell assay); PMID 10519380, Marutani et al. 1999 (TP53 yeast assay); PMID 11896595, Monti et al. 2002 (TP53 yeast assay). Items left out: 78.

| Gene | OS2 (functional evidence for) | SBS2 (against) | Clinical evidence, by best level |
|---|---|---|---|
| EGFR | 3 | 0 | A 2, B 2, C 11, D 12 |
| TP53 | 10 | 0 | B 5, C 5, D 8, E 1 |

## Classes

**EGFR**, by 3T3 score:

| 3T3 score | Oncogenic | Likely oncogenic | VUS | Likely benign | Benign | mixed | not assessed |
|---|---|---|---|---|---|---|---|
| 4 | 0 | 1 | 17 | 1 | 0 | 0 | 38 |
| 3 | 0 | 2 | 24 | 1 | 0 | 0 | 22 |
| 2 | 0 | 0 | 34 | 1 | 1 | 0 | 3 |
| 1 | 0 | 3 | 187 | 7 | 0 | 0 | 3 |

**TP53**, by Kato class:

| Kato class | Oncogenic | Likely oncogenic | VUS | Likely benign | Benign | mixed | not assessed |
|---|---|---|---|---|---|---|---|
| non-functional | 2 | 91 | 398 | 1 | 0 | 1 | 0 |
| partially functional | 0 | 5 | 571 | 5 | 0 | 1 | 0 |
| functional | 0 | 0 | 1123 | 21 | 0 | 3 | 0 |
| supertrans | 0 | 0 | 80 | 4 | 0 | 0 | 0 |

## Examples

| Gene | Variant | Assay | Class | Points | Criteria |
|---|---|---|---|---|---|
| EGFR | L858R | 3T3 4 | Likely oncogenic | 6 | OS3 +4, OP1 +1, OP4 +1 |
| EGFR | G719A | 3T3 3 | VUS | 4 | OM3 +2, OP1 +1, OP4 +1 |
| EGFR | G719C | 3T3 4 | VUS | 5 | OM4 +2, OP3 +1, OP1 +1, OP4 +1 |
| EGFR | G719S | 3T3 3 | VUS | 5 | OM4 +2, OP3 +1, OP1 +1, OP4 +1 |
| EGFR | L861Q | 3T3 3 | VUS | 4 | OM3 +2, OP1 +1, OP4 +1 |
| EGFR | S768I | 3T3 4 | VUS | 4 | OM3 +2, OP1 +1, OP4 +1 |
| EGFR | E709K | 3T3 4 | VUS | 4 | OM3 +2, OP1 +1, OP4 +1 |
| EGFR | T790M | 3T3 3 | Likely oncogenic | 9 | OS3 +4, OS2 +4, OP1 +1 |
| EGFR | L747P | 3T3 4 | VUS | 2 | OP3 +1, OP4 +1 |
| EGFR | E746_A750del | 3T3 3 | Likely oncogenic | 6 | OS3 +4, OM2 +2 |
| TP53 | R175H | non-functional | Likely oncogenic | 6 | OS3 +4, OP1 +1, OP4 +1 |
| TP53 | R248Q | non-functional | Likely oncogenic | 9 | OS3 +4, OS2 +4, OP1 +1 |
| TP53 | R248W | non-functional | Likely oncogenic | 6 | OS3 +4, OP1 +1, OP4 +1 |
| TP53 | R273H | non-functional | VUS | 5 | OS3 +4, OP1 +1 |
| TP53 | G245S | non-functional | Oncogenic | 10 | OS3 +4, OS2 +4, OP1 +1, OP4 +1 |
| TP53 | Y220C | non-functional | Likely oncogenic | 6 | OS3 +4, OP1 +1, OP4 +1 |
| TP53 | R337H | partially functional | VUS | 1 | OP3 +1 |
| TP53 | V217M | functional | VUS | 0 | OP4 +1, SBP1 -1 |
| TP53 | P47S | supertrans | Likely benign | -5 | SBS1 -4, SBP1 -1 |

Variants whose own allele has a ClinVar oncogenicity classification keep it in the ClinVar column; it never counts as evidence for itself. Insertions written as protein changes (most EGFR exon 20 insertions) cannot be traced to one DNA change and are not assessed.
