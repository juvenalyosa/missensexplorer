import io

import pandas as pd

from missensexplorer.inputs import parse_text, read_file, read_table, read_vcf, region_query


def kinds(text, gene=None):
    return [(v.kind, v.query) if v.ok else ("error", v.raw) for v in parse_text(text, gene=gene).variants]


def test_free_text_formats():
    got = kinds("""BRAF V600E
TP53: R175H, R248Q, R213*
chr7:140753336 A>T
7-140753336-A-T
7 140753336 A T
rs113488022
BRCA1:c.5096G>A
NM_004333.6(BRAF):c.1799T>A
NC_000007.14:g.140753336A>T
CFTR p.Phe508del
KRAS:p.(Gly12Asp)
""")
    assert got == [
        ("hgvs", "BRAF:p.Val600Glu"),
        ("hgvs", "TP53:p.Arg175His"),
        ("hgvs", "TP53:p.Arg248Gln"),
        ("hgvs", "TP53:p.Arg213Ter"),
        ("region", "7 140753336 . A T . . ."),
        ("region", "7 140753336 . A T . . ."),
        ("region", "7 140753336 . A T . . ."),
        ("id", "rs113488022"),
        ("hgvs", "BRCA1:c.5096G>A"),
        ("hgvs", "NM_004333.6:c.1799T>A"),
        ("hgvs", "NC_000007.14:g.140753336A>T"),
        ("hgvs", "CFTR:p.Phe508del"),
        ("hgvs", "KRAS:p.Gly12Asp"),
    ]


def test_default_gene_and_errors():
    assert kinds("V600E\nc.1799T>A", gene="BRAF") == [("hgvs", "BRAF:p.Val600Glu"), ("hgvs", "BRAF:c.1799T>A")]
    out = kinds("V600E\nthis is not a variant")
    assert out[0] == ("error", "V600E")
    assert all(k == "error" for k, _ in out)


def test_region_query_indels():
    assert region_query("7", 100, "A", "-") == "7 100 100 A/- 1"
    assert region_query("7", 100, "-", "T") == "7 101 100 -/T 1"
    assert region_query("7", 100, "ATCT", "A") == "7 100 . ATCT A . . ."


VCF = """##fileformat=VCFv4.2
##contig=<ID=chr1,length=248956422>
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	S1	S2
chr7	140753336	rs113488022	A	T	60	PASS	.	GT	0/1	0/0
chr17	7673802	.	C	T,A	60	PASS	.	GT	1/2	0/0
chr17	43063930	.	C	T	60	LowQ	.	GT	0/0	1/1
chr1	100	.	A	<DEL>	60	PASS	.	GT	0/1	0/1
chrX	5000	.	G	A	60	PASS	.	GT	1	./.
"""


def test_vcf_single_sample_zygosity_and_skips():
    p = read_vcf(io.StringIO(VCF), sample="S1")
    assert p.assembly_hint == "GRCh38"
    assert [(v.chrom, v.pos, v.alt, v.zygosity) for v in p.variants] == [
        ("7", 140753336, "T", "heterozygous"),
        ("17", 7673802, "T", "heterozygous"),
        ("17", 7673802, "A", "heterozygous"),
        ("X", 5000, "A", "hemizygous"),
    ]
    assert p.stats["skipped_not_carried"] == 1
    assert p.stats["skipped_symbolic"] == 1
    assert p.variants[0].vcf_id == "rs113488022"


def test_vcf_multi_sample_carriers_and_filters():
    p = read_vcf(io.StringIO(VCF), pass_only=True)
    assert all(v.carriers for v in p.variants)
    assert p.stats["skipped_filter"] == 1
    assert p.variants[-1].carriers == "S1:hemizygous"  # S2 has a missing call


def test_vcf_gene_regions_filter():
    p = read_vcf(io.StringIO(VCF), regions=[("17", 7661779, 7687546, "TP53")])
    assert {v.pos for v in p.variants} == {7673802}
    assert p.stats["skipped_outside_genes"] == 4


def test_read_file_dispatch_gz_and_table(tmp_path):
    import gzip

    f = tmp_path / "x.vcf.gz"
    f.write_bytes(gzip.compress(VCF.encode()))
    assert len(read_file(str(f), sample="S2").variants) == 1

    df = pd.DataFrame({"Chromosome": ["chr7"], "Position": [140753336], "Ref": ["A"], "Alt": ["T"]})
    assert read_table(df).variants[0].query == "7 140753336 . A T . . ."
    df2 = pd.DataFrame({"gene": ["TP53", "BRAF"], "mutation": ["R175H", "V600E"]})
    assert [v.query for v in read_table(df2).variants] == ["TP53:p.Arg175His", "BRAF:p.Val600Glu"]
