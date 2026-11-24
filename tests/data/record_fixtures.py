"""Re-record the VEP fixtures used by tests/test_edge_cases.py (needs internet).

python tests/data/record_fixtures.py
"""

import json
from pathlib import Path

from missensexplorer.api import ENSEMBL, Http
from missensexplorer.vep import BASE_PARAMS, PREDICTOR_PARAMS

HERE = Path(__file__).parent
PARAMS = {**BASE_PARAMS, **PREDICTOR_PARAMS}

CASES = {
    # name: (assembly, endpoint, payload key, inputs)
    "vep_edge_hgvs": (
        "GRCh38",
        "/vep/human/hgvs",
        "hgvs_notations",
        [
            "BRCA1:c.68_69del",  # frameshift (185delAG founder variant)
            "TP53:c.1A>G",  # start lost
            "TP53:c.375+1G>A",  # canonical splice donor
            "TP53:c.375G>A",  # synonymous, last base of exon 4 (splice-altering)
        ],
    ),
    "vep_edge_region": (
        "GRCh38",
        "/vep/human/region",
        "variants",
        [
            "MT 3460 . G A . . .",  # MT-ND1 A52T (LHON), mitochondrial codon table
            "X 31478302 . G A . . .",  # DMD, chromosome X
            "17 7673802 . C T . . .",  # TP53 R273H  } the two alleles of a multi-allelic VCF record
            "17 7673802 . C A . . .",  # TP53 R273L  }
        ],
    ),
    "vep_edge_grch37": (
        "GRCh37",
        "/vep/human/region",
        "variants",
        [
            "7 140453136 . A T . . .",  # BRAF V600E on GRCh37
        ],
    ),
}

# A protein change that two SNVs produce (CAA -> CAC or CAT), asked with ambiguous_hgvs: VEP answers once per
# transcript and allele, including alleles that change another residue on the reported transcript. Only the
# fields that vep.protein_matches reads are kept.
AMBIGUOUS = ("vep_ambiguous_kras_q61h", "KRAS:p.Gln61His")
KEEP = ("input", "seq_region_name", "start", "end", "allele_string", "strand")
KEEP_TC = ("gene_symbol", "transcript_id", "amino_acids", "protein_start", "consequence_terms")


def _trim(result):
    out = {k: result.get(k) for k in KEEP}
    out["transcript_consequences"] = [{k: tc.get(k) for k in KEEP_TC} for tc in result["transcript_consequences"]]
    return out


if __name__ == "__main__":
    http = Http()
    for name, (asm, path, key, inputs) in CASES.items():
        data = http.post_json(
            ENSEMBL[asm] + path, {key: inputs}, params=PARAMS, headers={"Content-Type": "application/json"}
        )
        (HERE / f"{name}.json").write_text(json.dumps(data, indent=1))
        print(name, len(data), "results")
    name, query = AMBIGUOUS
    data = http.post_json(
        ENSEMBL["GRCh38"] + "/vep/human/hgvs",
        {"hgvs_notations": [query]},
        params={**BASE_PARAMS, "ambiguous_hgvs": 1},
        headers={"Content-Type": "application/json"},
    )
    (HERE / f"{name}.json").write_text(json.dumps([_trim(r) for r in data], indent=1))
    print(name, len(data), "results")
