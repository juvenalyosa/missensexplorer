# Running VEP on your own computer

By default MissensExplorer uses the public Ensembl VEP REST API. The **local VEP** mode runs
[Ensembl VEP](https://www.ensembl.org/vep) on your computer instead. This is useful when:

- you analyse **exomes or genomes** (the public API is meant for thousands of variants, not millions);
- the data are **patient data** that must not leave the computer (combine it with *offline mode*);
- you need to work **without internet**.

> MissensExplorer talks to local VEP through standard input and output. That plumbing is covered by
> the test suite with a simulated VEP. Check your own installation with the test command in step 3.

## 1. Install Docker and the VEP image

Install [Docker Desktop](https://www.docker.com/products/docker-desktop/), then:

```bash
docker pull ensemblorg/ensembl-vep:release_116.0
mkdir -p $HOME/vep_data
```

## 2. Download the VEP cache and the genome (about 30 GB for GRCh38)

```bash
docker run --rm -it -v $HOME/vep_data:/data ensemblorg/ensembl-vep:release_116.0 \
  INSTALL.pl -a cf -s homo_sapiens -y GRCh38 -c /data
```

`-a cf` downloads the cache and the FASTA file. VEP needs the FASTA file to write HGVS notations.
Use `-y GRCh37` for hg19 data.

## 3. Check that VEP works

```bash
printf '7\t140753336\t140753336\tA/T\t+\ttest\n' | docker run --rm -i -v $HOME/vep_data:/data \
  ensemblorg/ensembl-vep:release_116.0 vep --cache --offline --dir_cache /data --format ensembl --json \
  --output_file STDOUT --hgvs --fasta /data/homo_sapiens/116_GRCh38/Homo_sapiens.GRCh38.dna.toplevel.fa.gz
```

The output should be one line of JSON that mentions `BRAF` and `p.Val600Glu`.

## 4. Tell MissensExplorer how to start VEP

Use the command from step 3, without the arguments that MissensExplorer adds itself (`--format`,
`--json`, `--output_file`, `--hgvs` and the other annotation flags).

**Desktop app:** open *Options → Annotation → Local VEP* and paste the command:

```
docker run --rm -i -v $HOME/vep_data:/data ensemblorg/ensembl-vep:release_116.0 vep --cache --offline --dir_cache /data --fasta /data/homo_sapiens/116_GRCh38/Homo_sapiens.GRCh38.dna.toplevel.fa.gz
```

**Command line:**

```bash
export MISSENSEXPLORER_VEP="docker run --rm -i -v $HOME/vep_data:/data ensemblorg/ensembl-vep:release_116.0 vep --cache --offline --dir_cache /data --fasta /data/homo_sapiens/116_GRCh38/Homo_sapiens.GRCh38.dna.toplevel.fa.gz"
missensexplorer -f patient.vcf.gz --vep local            # local VEP + online ClinVar/UniProt/AlphaFold
missensexplorer -f patient.vcf.gz --vep local --offline  # nothing leaves the computer
```

In local mode, VCF files and genomic coordinates (`chr7:140753336 A>T`) are annotated locally.
Variants typed as HGVS or rsIDs still go to the public API, because resolving them needs the
Ensembl database; in offline mode they are reported as *not analysed*.

## 5. Optional: predictors (plugins)

The local cache includes SIFT and PolyPhen-2. The other predictors are VEP plugins that read
their own data files. Download each file from its provider and check its licence (some are free
only for academic use). Then add the plugins to the extra arguments, adjusting the file names:

```bash
export MISSENSEXPLORER_VEP_ARGS="--dir_plugins /data/Plugins \
  --plugin AlphaMissense,file=/data/AlphaMissense_hg38.tsv.gz \
  --plugin CADD,snv=/data/whole_genome_SNVs.tsv.gz \
  --plugin dbNSFP,/data/dbNSFP_grch38.gz,REVEL_score \
  --plugin EVE,file=/data/eve_merged.vcf.gz \
  --plugin SpliceAI,snv=/data/spliceai_scores.raw.snv.hg38.vcf.gz,indel=/data/spliceai_scores.raw.indel.hg38.vcf.gz"
```

In the desktop app, paste the same text in *Extra arguments*. To install the plugin code, run
`INSTALL.pl -a p --PLUGINS AlphaMissense,CADD,dbNSFP,EVE,SpliceAI -c /data` with the same Docker
command as in step 2.

Without plugins, the verdict relies on SIFT and PolyPhen-2 alone, which is weaker evidence.
The *Run information* window lists which sources were used.
