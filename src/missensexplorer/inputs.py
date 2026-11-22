"""Turn whatever the user gives us into VEP-ready variant queries.

Accepted inputs
---------------
* VCF / VCF.gz files (single- or multi-sample, multi-allelic, SNVs and indels)
* CSV / TSV / Excel tables with chrom/pos/ref/alt columns, or a column of
  free-text variants (optionally with a gene column)
* Free text, one or many per line, separated by newlines, commas or semicolons:
    chr7:140753336 A>T     7-140753336-A-T     7 140753336 A T
    rs113488022
    BRAF V600E   BRAF:p.Val600Glu   TP53: R175H, R248Q, R273H
    BRCA1:c.5096G>A   NM_004333.6:c.1799T>A   NC_000007.14:g.140753336A>T
    V600E            (when a default gene is given)
"""

from __future__ import annotations

import gzip
import io
import re
from dataclasses import dataclass, field
from pathlib import Path

from .aa import AA1_TO_3, to_one_letter

CHROMS = {str(i) for i in range(1, 23)} | {"X", "Y", "MT"}


@dataclass
class VariantInput:
    raw: str                      # what the user wrote (or a file-derived description)
    kind: str = ""                # 'region' | 'hgvs' | 'id'
    query: str = ""               # exact string sent to VEP
    source: str = "text"
    chrom: str | None = None
    pos: int | None = None
    ref: str | None = None
    alt: str | None = None
    vcf_id: str | None = None
    vcf_filter: str | None = None
    sample: str | None = None
    genotype: str | None = None
    zygosity: str | None = None
    carriers: str | None = None
    error: str | None = None

    @property
    def ok(self) -> bool:
        return self.error is None and bool(self.query)


@dataclass
class ParsedInput:
    variants: list[VariantInput] = field(default_factory=list)
    samples: list[str] = field(default_factory=list)
    assembly_hint: str | None = None
    stats: dict = field(default_factory=dict)
    warnings: list[str] = field(default_factory=list)

    def extend(self, other: "ParsedInput") -> None:
        self.variants.extend(other.variants)
        for s in other.samples:
            if s not in self.samples:
                self.samples.append(s)
        self.assembly_hint = self.assembly_hint or other.assembly_hint
        for k, v in other.stats.items():
            self.stats[k] = self.stats.get(k, 0) + v
        self.warnings.extend(other.warnings)


# --------------------------------------------------------------------------- helpers

def norm_chrom(c: str) -> str:
    c = c.strip()
    if c.lower().startswith("chr"):
        c = c[3:]
    c = c.upper()
    if c == "M":
        c = "MT"
    return c


def region_query(chrom: str, pos: int, ref: str, alt: str) -> str:
    """Build a VEP region string. VCF-style when possible, Ensembl style for '-' alleles."""
    ref, alt = ref.upper(), alt.upper()
    if ref != "-" and alt != "-":
        return f"{chrom} {pos} . {ref} {alt} . . ."
    if alt == "-":  # deletion of `ref` starting at pos
        return f"{chrom} {pos} {pos + len(ref) - 1} {ref}/- 1"
    # insertion of `alt` after pos
    return f"{chrom} {pos + 1} {pos} -/{alt} 1"


_RS = re.compile(r"^rs\d+$", re.I)
_HGVS = re.compile(r"^(?P<ref>[A-Za-z0-9_.\-]+):(?P<type>[cgmnpr])\.(?P<rest>\S+)$")
_GENOMIC = re.compile(
    r"^(?:chr)?(?P<chrom>\d{1,2}|X|Y|M|MT)[\s:_\-]+(?P<pos>\d[\d,]*)[\s:_\-]*"
    r"(?P<ref>[ACGTN]+|-)[\s:_\-/>]+(?P<alt>[ACGTN]+|-)$",
    re.I,
)
_AA = r"(?:[A-Z][a-z]{2}|[ACDEFGHIKLMNPQRSTVWYUOX*])"
_PROT = re.compile(
    rf"^(?:p\.)?\(?(?P<ref>{_AA})(?P<pos>\d+)(?P<alt>{_AA}|=|Ter|del|dup|fs\S*|delins\S+|ins\S+|_\S+)\)?$"
)
_NUC_CHANGE = re.compile(r"^(?P<type>[cgmnr])\.(?P<rest>\S+)$")
_GENE = re.compile(r"^(?:[A-Za-z][A-Za-z0-9\-]{0,19}(?:orf\d+)?|ENS[TPG]\d+(?:\.\d+)?|N[MPRC]_\d+(?:\.\d+)?)$")


def _protein_change(change: str) -> str | None:
    """'V600E' / 'p.Val600Glu' / 'R213*' / 'R213X' -> 'p.Val600Glu' (VEP-friendly), else None."""
    m = _PROT.match(change.strip())
    if not m:
        return None
    ref, pos, alt = m.group("ref"), m.group("pos"), m.group("alt")
    r1 = to_one_letter(ref)
    if r1 is None or r1 == "*":
        return None
    ref3 = AA1_TO_3[r1]
    if alt in ("X", "*", "Ter"):
        alt3 = "Ter"
    elif alt == "=" or alt.startswith(("del", "dup", "fs", "ins", "_")):
        alt3 = alt
    else:
        a1 = to_one_letter(alt)
        if a1 is None:
            return None
        alt3 = AA1_TO_3[a1]
    return f"p.{ref3}{pos}{alt3}"


def _clean_hgvs(s: str) -> str:
    s = re.sub(r"\(([A-Za-z0-9\-]+)\)(?=:)", "", s)          # NM_004333.6(BRAF):c. -> NM_004333.6:c.
    s = re.sub(r"^chr", "", s, flags=re.I)                      # chr7:g. -> 7:g.
    s = re.sub(r":p\.\((.+)\)$", r":p.\1", s)                  # :p.(Val600Glu) -> :p.Val600Glu
    return s


def parse_token(token: str, gene: str | None = None) -> VariantInput | None:
    """Parse one variant expression. Returns None if it is not a variant (e.g. a bare gene)."""
    t = token.strip().strip(",;")
    if not t:
        return None
    if _RS.match(t):
        return VariantInput(raw=t, kind="id", query=t.lower())

    m = _GENOMIC.match(t)
    if m:
        chrom = norm_chrom(m.group("chrom"))
        pos = int(m.group("pos").replace(",", ""))
        ref, alt = m.group("ref").upper(), m.group("alt").upper()
        return VariantInput(raw=t, kind="region", query=region_query(chrom, pos, ref, alt),
                            chrom=chrom, pos=pos, ref=ref, alt=alt)

    h = _clean_hgvs(t)
    m = _HGVS.match(h)
    if m:
        if m.group("type") == "p":
            pc = _protein_change(m.group("rest"))
            if pc:
                h = f"{m.group('ref')}:{pc}"
        return VariantInput(raw=t, kind="hgvs", query=h)

    # "GENE CHANGE" / "GENE:CHANGE" written in one token pair
    gm = re.match(r"^(?P<gene>[A-Za-z0-9_.\-]+)[\s:]+(?P<change>\S+)$", t)
    if gm and _GENE.match(gm.group("gene")):
        sub = parse_token(gm.group("change"), gene=gm.group("gene"))
        if sub is not None and sub.kind == "hgvs":
            sub.raw = t
            return sub

    pc = _protein_change(t)
    nc = _NUC_CHANGE.match(t)
    if pc or nc:
        if not gene:
            return VariantInput(raw=t, error="needs a gene, e.g. 'BRAF V600E' (or set the default gene)")
        return VariantInput(raw=t if gene in t else f"{gene} {t}", kind="hgvs",
                            query=f"{gene}:{pc}" if pc else f"{gene}:{t}")
    return None


def _looks_like_vcf_line(fields: list[str]) -> bool:
    return len(fields) >= 5 and fields[1].isdigit() and re.fullmatch(r"[ACGTNacgtn]+", fields[3] or "") is not None


def parse_text(text: str, gene: str | None = None, source: str = "text") -> ParsedInput:
    """Parse free text: one or more variants per line."""
    out = ParsedInput()
    gene = gene.strip() if gene and gene.strip() else None
    for raw_line in text.splitlines():
        line = raw_line.split("#", 1)[0].strip() if not raw_line.lstrip().startswith("#") else ""
        if not line:
            continue
        fields = line.split("\t")
        if _looks_like_vcf_line(fields):  # a pasted VCF data line
            chrom = norm_chrom(fields[0])
            for alt in fields[4].split(","):
                if re.fullmatch(r"[ACGTNacgtn]+", alt):
                    v = VariantInput(raw=f"{fields[0]}:{fields[1]} {fields[3]}>{alt}", kind="region",
                                     query=region_query(chrom, int(fields[1]), fields[3], alt),
                                     chrom=chrom, pos=int(fields[1]), ref=fields[3].upper(), alt=alt.upper(),
                                     vcf_id=fields[2] if fields[2] != "." else None, source=source)
                    out.variants.append(v)
            continue

        # "GENE: change, change ..." -> the gene carries over to the following items
        current_gene = gene
        for piece in re.split(r"[,;]", line):
            piece = piece.strip()
            if not piece:
                continue
            v = parse_token(piece, current_gene)
            found = [v] if v is not None else []
            if v is None:
                for tok in piece.split():
                    tv = parse_token(tok, current_gene)
                    if tv is None:
                        cand = tok.rstrip(":")
                        if _GENE.match(cand):
                            current_gene = cand
                            continue
                        tv = VariantInput(raw=tok, error="unrecognised variant format")
                    found.append(tv)
            for fv in found:
                if fv.kind == "hgvs" and fv.ok:
                    prefix = fv.query.split(":", 1)[0]
                    if _GENE.match(prefix) and not prefix.startswith(("NC_", "NM_", "NP_", "NR_", "ENS")):
                        current_gene = prefix
                fv.source = source
                out.variants.append(fv)
    out.stats["parsed"] = sum(v.ok for v in out.variants)
    out.stats["unparsed"] = sum(not v.ok for v in out.variants)
    return out


# --------------------------------------------------------------------------- VCF

_GRCH37_LEN = {"1": 249250621, "2": 243199373, "X": 155270560}
_GRCH38_LEN = {"1": 248956422, "2": 242193529, "X": 156040895}


def _assembly_from_header(line: str) -> str | None:
    low = line.lower()
    if line.startswith("##contig"):
        m_id = re.search(r"ID=([^,>]+)", line)
        m_len = re.search(r"length=(\d+)", line)
        if m_id and m_len:
            c, n = norm_chrom(m_id.group(1)), int(m_len.group(1))
            if _GRCH38_LEN.get(c) == n:
                return "GRCh38"
            if _GRCH37_LEN.get(c) == n:
                return "GRCh37"
    if line.startswith(("##reference", "##contig", "##assembly")):
        if re.search(r"grch38|hg38|hs38|b38", low):
            return "GRCh38"
        if re.search(r"grch37|hg19|hs37|b37|human_g1k_v37", low):
            return "GRCh37"
    return None


def _zygosity(alleles: list[str], idx: str) -> str | None:
    called = [a for a in alleles if a not in (".", "")]
    if not called:  # no-call: the sample is not known to carry the allele
        return None
    n = called.count(idx)
    if n == 0:
        return None
    if len(alleles) == 1:
        return "hemizygous"
    if n == len(alleles):
        return "homozygous"
    return "heterozygous"


def read_vcf(lines, source: str = "vcf", sample: str | None = None,
             regions: list[tuple[str, int, int, str]] | None = None,
             pass_only: bool = False) -> ParsedInput:
    """Stream a VCF. `regions` = [(chrom, start, end, gene)] keeps only variants inside them."""
    out = ParsedInput()
    st = {"vcf_records": 0, "skipped_not_carried": 0, "skipped_symbolic": 0,
          "skipped_outside_genes": 0, "skipped_filter": 0}
    samples: list[str] = []
    sample_idx: list[int] = []
    reg_by_chrom: dict[str, list[tuple[int, int]]] = {}
    for c, s, e, _ in regions or []:
        reg_by_chrom.setdefault(norm_chrom(c), []).append((s, e))

    for line in lines:
        if isinstance(line, bytes):
            line = line.decode("utf-8", "replace")
        line = line.rstrip("\r\n")
        if not line:
            continue
        if line.startswith("##"):
            out.assembly_hint = out.assembly_hint or _assembly_from_header(line)
            continue
        if line.startswith("#"):
            hdr = line.lstrip("#").split("\t")
            samples = hdr[9:] if len(hdr) > 9 else []
            out.samples = samples
            if sample:
                if sample not in samples:
                    raise ValueError(f"sample '{sample}' not in VCF (available: {', '.join(samples) or 'none'})")
                sample_idx = [samples.index(sample)]
            elif len(samples) == 1:
                sample_idx = [0]
            else:
                sample_idx = list(range(len(samples)))
            continue

        f = line.split("\t")
        if len(f) < 5:
            continue
        st["vcf_records"] += 1
        chrom = norm_chrom(f[0])
        try:
            pos = int(f[1])
        except ValueError:
            continue
        if regions is not None:
            regs = reg_by_chrom.get(chrom)
            if not regs or not any(s <= pos <= e for s, e in regs):
                st["skipped_outside_genes"] += 1
                continue
        flt = f[6] if len(f) > 6 else "."
        if pass_only and flt not in ("PASS", "."):
            st["skipped_filter"] += 1
            continue
        ref = f[3].upper()
        alts = f[4].split(",")
        fmt = f[8].split(":") if len(f) > 8 else []
        gt_i = fmt.index("GT") if "GT" in fmt else None
        gts = []
        if gt_i is not None:
            for i in sample_idx:
                col = f[9 + i].split(":") if 9 + i < len(f) else []
                gts.append((samples[i], col[gt_i] if gt_i < len(col) else "."))

        for ai, alt in enumerate(alts, start=1):
            alt = alt.upper()
            if alt in (".", "*") or alt.startswith("<") or "[" in alt or "]" in alt or not re.fullmatch(r"[ACGTN]+", alt) or not re.fullmatch(r"[ACGTN]+", ref):
                st["skipped_symbolic"] += 1
                continue
            v = VariantInput(raw=f"{f[0]}:{pos} {ref}>{alt}", kind="region",
                             query=region_query(chrom, pos, ref, alt), source=source,
                             chrom=chrom, pos=pos, ref=ref, alt=alt,
                             vcf_id=f[2] if f[2] not in (".", "") else None, vcf_filter=flt)
            if gts:
                carriers = []
                for sname, gt in gts:
                    z = _zygosity(re.split(r"[/|]", gt), str(ai))
                    if z:
                        carriers.append((sname, gt, z))
                if not carriers:
                    st["skipped_not_carried"] += 1
                    continue
                if len(sample_idx) == 1:
                    v.sample, v.genotype, v.zygosity = carriers[0]
                else:
                    v.carriers = "; ".join(f"{s}:{z}" for s, _, z in carriers)
            out.variants.append(v)

    st["variants"] = len(out.variants)
    out.stats = st
    return out


# --------------------------------------------------------------------------- tables

_COLS = {
    "chrom": ("chrom", "#chrom", "chr", "chromosome", "chro", "seqname", "contig"),
    "pos": ("pos", "position", "start", "bp", "genomic position", "genomic_position"),
    "ref": ("ref", "reference", "ref_allele", "ref_alle", "reference_allele"),
    "alt": ("alt", "alternate", "alt_allele", "alt_alle", "alt_alleles", "alternative", "allele"),
    "text": ("variant", "variants", "hgvs", "mutation", "mutations", "input", "rsid", "rs", "id",
             "change", "protein_change", "hgvsp", "hgvsc"),
    "gene": ("gene", "symbol", "gene_symbol", "hugo_symbol"),
}


def read_table(df, source: str = "table", gene: str | None = None) -> ParsedInput:
    cols = {c.lower().strip(): c for c in df.columns}

    def pick(key):
        for alias in _COLS[key]:
            if alias in cols:
                return cols[alias]
        return None

    c_chr, c_pos, c_ref, c_alt = pick("chrom"), pick("pos"), pick("ref"), pick("alt")
    out = ParsedInput()
    if c_chr and c_pos and c_ref and c_alt:
        for _, row in df.iterrows():
            try:
                chrom, pos = norm_chrom(str(row[c_chr])), int(row[c_pos])
            except (TypeError, ValueError):
                continue
            ref, alt = str(row[c_ref]).strip().upper(), str(row[c_alt]).strip().upper()
            for a in alt.split(","):
                if re.fullmatch(r"[ACGTN]+|-", ref) and re.fullmatch(r"[ACGTN]+|-", a):
                    out.variants.append(VariantInput(raw=f"{row[c_chr]}:{pos} {ref}>{a}", kind="region",
                                                     query=region_query(chrom, pos, ref, a), source=source,
                                                     chrom=chrom, pos=pos, ref=ref, alt=a))
        return out
    c_txt, c_gene = pick("text"), pick("gene")
    if c_txt is None:
        c_txt = df.columns[0]
    for _, row in df.iterrows():
        val = row[c_txt]
        if val is None or (isinstance(val, float) and val != val):
            continue
        g = str(row[c_gene]).strip() if c_gene and isinstance(row[c_gene], str) else gene
        out.extend(parse_text(str(val), gene=g, source=source))
    return out


# --------------------------------------------------------------------------- dispatcher

def read_file(name: str, data: bytes | None = None, *, gene: str | None = None,
              sample: str | None = None, regions=None, pass_only: bool = False) -> ParsedInput:
    """Read a VCF / VCF.gz / CSV / TSV / XLSX / TXT file given its name and (optionally) bytes.

    VCF files given by path are streamed line by line (never fully loaded in memory).
    """
    if data is None:
        with open(name, "rb") as fh:
            head = fh.read(2)
            gz = head == b"\x1f\x8b"
        opener = gzip.open if gz else open
        with opener(name, "rt", encoding="utf-8", errors="replace") as fh:
            first = fh.read(4096)
        low = Path(name).name.lower().removesuffix(".gz").removesuffix(".bgz")
        if low.endswith(".vcf") or first.startswith("##fileformat=VCF") or first.startswith("#CHROM"):
            with opener(name, "rt", encoding="utf-8", errors="replace") as fh:
                return read_vcf(fh, source=Path(name).name, sample=sample, regions=regions, pass_only=pass_only)
        data = Path(name).read_bytes()
    base = Path(name).name
    lower = base.lower()
    if data[:2] == b"\x1f\x8b":
        data = gzip.decompress(data)
        lower = lower.removesuffix(".gz").removesuffix(".bgz")
    head = data[:4096].decode("utf-8", "replace")
    if lower.endswith(".vcf") or head.startswith("##fileformat=VCF") or head.startswith("#CHROM"):
        return read_vcf(io.StringIO(data.decode("utf-8", "replace")), source=base, sample=sample,
                        regions=regions, pass_only=pass_only)
    if lower.endswith((".xlsx", ".xls")):
        import pandas as pd
        return read_table(pd.read_excel(io.BytesIO(data)), source=base, gene=gene)
    if lower.endswith((".csv", ".tsv", ".tab")):
        import pandas as pd
        sep = "," if lower.endswith(".csv") else "\t"
        df = pd.read_csv(io.BytesIO(data), sep=sep, dtype=str, comment=None)
        if len(df.columns) > 1:
            return read_table(df, source=base, gene=gene)
    return parse_text(data.decode("utf-8", "replace"), gene=gene, source=base)


def vcf_samples(name: str, data: bytes) -> list[str]:
    """Sample names from a VCF header (reads only the header)."""
    if data[:2] == b"\x1f\x8b":
        import zlib
        d = zlib.decompressobj(16 + zlib.MAX_WBITS)
        data = d.decompress(data[:4_000_000], 8_000_000)
    for line in data.decode("utf-8", "replace").splitlines():
        if line.startswith("#CHROM"):
            return line.split("\t")[9:]
        if not line.startswith("#"):
            break
    return []
