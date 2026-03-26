# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Overview

This is a general-purpose Snakemake workflow for RNA-seq processing and differential expression analysis. It supports both **paired-end** and **single-end** reads, takes raw FASTQ files, and produces pairwise DEG (Differentially Expressed Genes) results that can be used for any downstream analysis — including as direct input for libncker's `exclusive` and `neighbors` modules. Designed to run on an SGE (Sun Grid Engine) HPC cluster with environment modules, but also works locally.

## Local Development Setup

The workflow uses `module load` on the HPC cluster. For local testing, create a venv once:

```bash
python3 -m venv snakemake/.venv
source snakemake/.venv/bin/activate
pip install pandas numpy matplotlib pydeseq2
```

Activate it before running any script locally: `source snakemake/.venv/bin/activate`

## Running the Workflow

```bash
# Dry run (check what will be executed)
snakemake -n

# Run locally
snakemake --cores <N>

# Submit to SGE cluster
snakemake --cluster "qsub {cluster.options}" --cluster-config cluster.json --jobs <N>

# Run a specific rule
snakemake <output_file>

# Run DESeq2 — included in rule all, runs automatically
snakemake --cores <N>

# Clean all outputs
snakemake clean
```

## Configuration (`config.yaml`)

| Key | Purpose |
|---|---|
| `aligner` | Aligner to use: `hisat2`, `bowtie2`, or `bwa` |
| `layout` | Read layout: `paired` (default) or `single` |
| `paths.data_dir` | Directory containing raw FASTQs — paired: `{sample}_1.fastq.gz`/`{sample}_2.fastq.gz`; single: `{sample}.fastq.gz` |
| `paths.output_dir` | Root output directory — all subdirs are created automatically |
| `paths.genome` | Path to reference genome FASTA |
| `paths.annotation` | Optional path to annotation GTF; leave empty if not available |
| `paths.splice_sites` | Optional pre-generated splice sites file (hisat2 only) |
| `python_bin` | Python interpreter to use for scripts — set to `.venv/bin/python3` on workstation, `python3` on HPC |
| `modules.bowtie2` / `modules.bwa` | Module names on the cluster — ignored on workstations |
| `modules.python` | Python module — must have `pandas`, `pydeseq2`, `matplotlib` installed |
| `metadata` | Path to `metadata.csv` for DESeq2 (see format below) |
| `deseq2.padj` | Adjusted p-value cutoff for significance (default: `0.05`) |
| `deseq2.log2fc` | Minimum absolute log2 fold change (default: `1.5`) |
| `deseq2.min_counts` | Minimum total counts across samples to retain a gene (default: `10`) |

## Pipeline Architecture

The workflow processes FASTQ files through these stages (in dependency order):

1. **check_pairs** *(paired only)* — `seqkit stats` to detect R1/R2 count mismatches
2. **repair_files** *(paired only)* — if counts differ, runs `bbmap repair.sh`; otherwise creates symlinks
3. **check_reads** *(single only)* — `seqkit stats` on the single read file
4. **fastqc** — per-sample QC; paired runs on R1+R2, single on the one file
5. **multiqc** — aggregates all FastQC ZIP reports
6. **build_index** — builds aligner index on the reference genome (one-time)
7. **align** — splice-aware or standard alignment; outputs a temp SAM + alignment summary
8. **sam_to_bam** — `samtools view` (temp BAM)
9. **sort_bam** — `samtools sort` → `{sample}.sorted.bam`
10. **flagstat** — `samtools flagstat` on sorted BAM
11. **coverage** — `samtools coverage` per-contig: extracts columns `rname` + `numreads`
12. **merge_counts** — `scripts/merge_counts.py` merges all `.coverage` files into a genes × samples count table
13. **deseq2** — `scripts/deseq2.py` runs all pairwise DESeq2 comparisons; produces one `<condAvscondB>_intersect.txt` per pair, full results, normalized counts, and MA plots

Intermediate SAM and unsorted BAM files are marked `temp()` and deleted automatically after downstream rules complete.

### Layout differences
- **Paired-end:** reads go through `repair_files` (sync/repair step) before alignment; aligned with `-1`/`-2` flags
- **Single-end:** no repair step; reads aligned directly from `data_dir` with `-U` flag (hisat2/bowtie2) or as a single positional arg (bwa); hisat2 strandness changes from `RF` to `R`

### Aligner notes
- **hisat2** — splice-aware (recommended for RNA-seq); also outputs novel splice sites
- **bowtie2** — not splice-aware; valid for organisms without known splice sites
- **bwa** — not splice-aware; writes BWA index alongside the genome file (`.amb`, `.ann`, `.bwt`, `.pac`, `.sa`). Separate indexes are built for hisat2 and bowtie2 (different output dirs).

All output subdirectories (`fastqc_reports/`, `processed_bams/`, `stats/`, `coverage/`, `counts/`, `deseq2/`) are created automatically under `paths.output_dir`. Index prefixes are also placed there.

Sample names are auto-detected from filenames in `paths.data_dir`: paired matches `{sample}_1.fastq.gz` / `{sample}_2.fastq.gz`; single matches `{sample}.fastq.gz` (files ending in `_1.fastq.gz` or `_2.fastq.gz` are excluded from single detection).

## Cluster Configuration (`cluster.json`)

SGE options per rule. The `__default__` applies to rules not explicitly listed:

| Rule | Threads | Memory |
|---|---|---|
| default | 1 | 10G |
| check_pairs, build_index | 6 | 16G |
| align, sam_to_bam, sort_bam | 8 | 32G / 16G / 10G |
| flagstat | 4 | 10G |
| merge_counts | 1 | 8G |
| deseq2 | 1 | 32G |

Thread count is accessed inside shell blocks via `$NSLOTS` (SGE environment variable).

## Environment Modules

Each shell rule loads its own modules via `module load`. The pattern used throughout:
```bash
set +u
source ~/.bashrc
module load programs/<tool-version>
```

Each shell block uses `command -v module &>/dev/null && module load <tool>` — on HPC the module is loaded; on a workstation the line is silently skipped and the tool is expected on `$PATH`. Thread count uses `${NSLOTS:-$(nproc)}` — `$NSLOTS` on SGE, `nproc` on workstation.

Modules used on HPC: `seqkit-2.5.1`, `bbmap-38.84`, `hisat2-2.1.0`, `fastqc-0.12.1`, `multiqc-1.23`, `samtools-1.10`. Bowtie2, BWA, and Python module names are configured in `config.yaml`.

## Scripts

| Script | Language | Purpose |
|---|---|---|
| `scripts/merge_counts.py` | Python | Reads all `{sample}.coverage` files, merges into a tab-separated genes × samples count table. Skips the `#rname` header line, fills missing values with 0, drops all-zero rows. |
| `scripts/deseq2.py` | Python | Detects all conditions from metadata, runs PyDESeq2 for every pairwise combination. Per pair outputs: `<condAvscondB>_intersect.txt` (ID + Expression, sorted by ID — libncker input), full results, normalized counts, MA plot. One shared PCA plot (all samples, log1p raw counts). Thresholds (padj, log2fc, min_counts) passed as CLI args from Snakemake params. |

## DESeq2 Metadata Format

Edit `metadata.csv` before running the `deseq2` rule. Column `sample` must match column names in `count_table.txt`:

```
sample,condition
SRR001,control
SRR002,control
SRR003,treatment
SRR004,treatment
```

## Notes

- The `splice_sites` rule is commented out — splice site extraction from the GTF was done separately; the file is pre-generated at `splice_sites_file`.
- HISAT2 alignment uses `--rna-strandness RF --dta` (reverse-stranded, downstream transcriptome assembly mode).
- The `annotation_gtf` path is checked at workflow load time (`USE_GTF = os.path.isfile(annotation_gtf)`), but this flag is not currently wired into any active rule.
- The `align` rule replaced the old `hisat2` rule — update any saved SGE job names accordingly.
