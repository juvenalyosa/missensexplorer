"""Run Ensembl VEP on this computer instead of the public REST API.

Why: no request limits, whole genomes are feasible, works without internet once the VEP
cache is installed, and patient variants never leave the machine.

The command is configurable, for example:

    vep                                       (VEP installed with its cache)
    docker run --rm -i -v $HOME/vep_data:/data ensemblorg/ensembl-vep vep --dir_cache /data

Variants are sent on standard input in VEP's default format and results are read as JSON
from standard output, so the same command works for a native install and for Docker.
Plugin options (AlphaMissense, CADD, dbNSFP/REVEL, SpliceAI, EVE need local data files) go
in the extra arguments, e.g. ``--plugin AlphaMissense,file=/data/AlphaMissense_hg38.tsv.gz``.
"""

from __future__ import annotations

import hashlib
import json
import logging
import os
import re
import shlex
import subprocess

log = logging.getLogger(__name__)

# Flags matching what the REST client requests (BASE_PARAMS) plus what REST adds by default.
LOCAL_FLAGS = [
    "--format", "ensembl", "--json", "--no_stats", "--force_overwrite",
    "--pick_allele_gene", "--hgvs", "--mane", "--canonical", "--uniprot", "--numbers", "--variant_class",
    "--protein", "--symbol", "--sift", "b", "--polyphen", "b", "--check_existing", "--af_gnomade",
    "--af_gnomadg", "--pubmed", "--var_synonyms",
]

ENV_COMMAND = "MISSENSEXPLORER_VEP"
ENV_ARGS = "MISSENSEXPLORER_VEP_ARGS"


def default_command() -> str:
    return os.environ.get(ENV_COMMAND, "vep --cache --offline")


def default_args() -> str:
    return os.environ.get(ENV_ARGS, "")


def to_default_format(query: str, ident: str) -> str | None:
    """Region query (VCF-like '7 140753336 . A T . . .' or Ensembl '7 100 100 A/- 1')
    -> VEP default format line 'chrom start end ref/alt strand id'."""
    parts = query.split()
    if len(parts) >= 5 and "/" in parts[3]:  # already Ensembl default format
        chrom, start, end, allele = parts[0], parts[1], parts[2], parts[3]
        strand = parts[4] if len(parts) > 4 else "1"
        return f"{chrom}\t{start}\t{end}\t{allele}\t{strand}\t{ident}"
    if len(parts) < 5:
        return None
    chrom, ref, alt = parts[0], parts[3].upper(), parts[4].upper()
    try:
        pos = int(parts[1])
    except ValueError:
        return None
    # trim the shared anchor base(s), as VEP does for VCF input
    while ref and alt and ref[0] == alt[0] and (len(ref) > 1 or len(alt) > 1):
        ref, alt, pos = ref[1:], alt[1:], pos + 1
    while ref and alt and ref[-1] == alt[-1] and (len(ref) > 1 or len(alt) > 1):
        ref, alt = ref[:-1], alt[:-1]
    if not ref:      # insertion between pos-1 and pos
        return f"{chrom}\t{pos}\t{pos - 1}\t-/{alt}\t+\t{ident}"
    end = pos + len(ref) - 1
    return f"{chrom}\t{pos}\t{end}\t{ref}/{alt or '-'}\t+\t{ident}"


class LocalVep:
    name = "Ensembl VEP (local)"

    def __init__(self, cache, assembly: str = "GRCh38", command: str | None = None, extra_args: str | None = None,
                 status=None, timeout: float | None = None, chunk: int = 5000):
        self.cache, self.assembly, self.status = cache, assembly, status
        self.command = command or default_command()
        self.extra_args = extra_args if extra_args is not None else default_args()
        self.timeout = timeout
        self.chunk = chunk
        self.release = None
        self._sig = hashlib.sha1(f"{self.command}|{self.extra_args}|{assembly}".encode()).hexdigest()[:10]

    def argv(self) -> list[str]:
        return (shlex.split(os.path.expandvars(os.path.expanduser(self.command))) + LOCAL_FLAGS
                + ["--assembly", self.assembly, "--output_file", "STDOUT"]
                + shlex.split(os.path.expandvars(self.extra_args)))

    def version(self) -> str | None:
        """'ensembl-vep 115.0' from `vep --help`, or None."""
        try:
            p = subprocess.run(shlex.split(os.path.expanduser(self.command)) + ["--help"], capture_output=True,
                               text=True, timeout=120)
        except (OSError, subprocess.TimeoutExpired):
            return None
        m = re.search(r"ensembl-vep\s*:\s*([\d.]+)", p.stdout + p.stderr)
        return f"ensembl-vep {m.group(1)}" if m else None

    def _key(self, q: str) -> str:
        return f"vep-local|{self._sig}|{q}"

    def annotate(self, queries: dict[str, list[str]], params=None, progress=None, stage: str = "VEP (local)",
                 source: str = "Ensembl VEP (local)") -> dict[tuple[str, str], dict]:
        out: dict[tuple[str, str], dict] = {}
        for kind, qs in queries.items():
            if kind != "region":
                for q in qs:
                    out[(kind, q)] = {"error": "local VEP needs genomic coordinates (VCF or chrom:pos ref>alt); "
                                               "HGVS and rsIDs need the online mode", "cache_error": True}
        qs = list(dict.fromkeys(queries.get("region", [])))
        if not qs:
            return out
        if self.status is not None:
            self.status.use(source, len(qs))
        hit = self.cache.get_many([self._key(q) for q in qs])
        todo = [q for q in qs if self._key(q) not in hit]
        for q in qs:
            if self._key(q) in hit:
                out[("region", q)] = hit[self._key(q)]
        total, done = len(todo), 0
        if progress:
            progress(stage, 0, total)
        for i in range(0, len(todo), self.chunk):
            batch = todo[i: i + self.chunk]
            res = self._run(batch, source)
            cacheable = {}
            for q, val in res.items():
                out[("region", q)] = val
                if "results" in val or val.get("cache_error"):
                    cacheable[self._key(q)] = val
            self.cache.set_many(cacheable)
            done += len(batch)
            if progress:
                progress(stage, done, total)
        return out

    def _run(self, batch: list[str], source: str) -> dict[str, dict]:
        lines, ids = [], {}
        res: dict[str, dict] = {}
        for n, q in enumerate(batch):
            ident = f"mx{n}"
            line = to_default_format(q, ident)
            if line is None:
                res[q] = {"error": "could not convert to VEP input", "cache_error": True}
                continue
            ids[ident] = q
            lines.append(line)
        if not lines:
            return res
        argv = self.argv()
        try:
            p = subprocess.run(argv, input="\n".join(lines) + "\n", capture_output=True, text=True,
                               timeout=self.timeout)
        except (OSError, subprocess.TimeoutExpired) as e:
            msg = f"could not run local VEP ({argv[0]}): {e}"
            return self._fail(batch, res, msg, source)
        if p.returncode != 0:
            tail = " ".join((p.stderr or p.stdout).strip().splitlines()[-3:])[:400]
            return self._fail(batch, res, f"local VEP exited with code {p.returncode}: {tail}", source)
        grouped: dict[str, list] = {}
        for line in p.stdout.splitlines():
            line = line.strip()
            if not line.startswith("{"):
                continue
            try:
                item = json.loads(line)
            except ValueError:
                continue
            q = ids.get(str(item.get("id")))
            if q is not None:
                grouped.setdefault(q, []).append(item)
        for ident, q in ids.items():
            res[q] = ({"results": grouped[q]} if q in grouped else
                      {"error": "VEP could not map this variant (check the assembly and the alleles)", "cache_error": True})
        return res

    def _fail(self, batch, res, msg, source):
        log.warning(msg)
        if self.status is not None:
            self.status.fail(source, msg, batch)
        for q in batch:
            res.setdefault(q, {"error": msg, "service_error": True})
        return res
