# libncker

Libncker is a lightweight CLI to:
1) identify tissue-exclusive lncRNAs from differential expression intersect tables,
2) compute a cis-neighborhood report (K nearest mRNAs) plus a simple `cis_score`,
3) summarize intergenic or intron-length distributions with geometric means, and
4) annotate differentially expressed mRNAs with NCBI gene summaries.

---

## Table of contents
- [Install](#install)
- [Quickstart](#quickstart)
- [Commands](#commands)
  - [1) extract-ids](#1-extract-lncrna-ids-from-gff)
  - [2) exclusive](#2-compute-tissue-exclusive-lncrnas-from-intersect-tables)
  - [3) neighbors](#3-compute-neighbors--cis-module-output)
  - [4) run](#4-run-everything-end-to-end)
  - [5) intergenic-stats](#5-summarize-intergenic-distances-from-a-gff)
  - [6) intron-stats](#6-summarize-intron-lengths-from-a-tsv-column)
  - [7) annotate](#7-annotate-de-mrnas-with-ncbi-gene-summaries)
- [Inputs](#inputs)
- [Outputs](#outputs)
- [How scoring works](#how-scoring-works)
- [Common pitfalls](#common-pitfalls)
- [Contributing](#contributing)
- [License](#license)

---

## Install

**Installation note:**
On some Debian/Ubuntu systems, Python is marked as externally managed (PEP 668).
In that case, please install using a virtual environment or ensure python3-full is available.

### Option A (recommended, HPC-friendly): virtual environment
```bash
python3 -m venv .venv
source .venv/bin/activate
python3 -m pip install -U pip
python3 -m pip install "git+https://github.com/yahelgb/libncker.git"
libncker -h
```
### Option B: local clone
```bash
git clone https://github.com/yahelgb/libncker.git
cd libncker
python3 -m pip install .
libncker -h
```
### Option C: user install (no venv)
```bash
export PATH="$HOME/.local/bin:$PATH"
python3 -m pip install --user "git+https://github.com/yahelgb/libncker.git"
which libncker
libncker -h
```
### Option D: run without PATH
```bash
python3 -m pip install --user "git+https://github.com/yahelgb/libncker.git"
python3 -m libncker -h
```
---

## Quickstart

### Full pipeline (recommended: transcript IDs)
```bash
libncker run \
  --gff /path/to/genome.gff \
  --lnc-ids-from-gff \
  --level transcript \
  --intersects /path/to/*_intersect.txt \
  --outdir results \
  --mode strict \
  --k 5
```

### Full pipeline (gene IDs)
```bash
libncker run \
  --gff /path/to/genome.gff \
  --lnc-ids-from-gff \
  --level gene \
  --intersects /path/to/*_intersect.txt \
  --outdir results \
  --mode strict \
  --k 5
```

Compressed `.gff.gz` inputs are supported directly.
---

## Commands

### 1) Extract lncRNA IDs from GFF
```bash
libncker extract-ids \
  --gff /path/to/genome.gff \
  --level transcript \
  --out lncRNA_IDs.txt \
  --summary lncRNA_IDs.summary.tsv
```
### 2) Compute tissue-exclusive lncRNAs from intersect tables
```bash
libncker exclusive \
  --lnc-ids lncRNA_IDs.txt \
  --level transcript \
  --intersects /path/to/*_intersect.txt \
  --outdir results \
  --mode strict
```
Lenient mode (for incomplete pairwise sets or large N):
```bash
libncker exclusive \
  --lnc-ids lncRNA_IDs.txt \
  --level transcript \
  --intersects /path/to/*_intersect.txt \
  --outdir results \
  --mode lenient \
  --min-required 3
```
### 3) Compute neighbors + cis module output
```bash
libncker neighbors \
  --gff /path/to/genome.gff \
  --level transcript \
  --exclusive results/*_exclusive_lncRNAs_strict.txt \
  --intersects /path/to/*vs*_intersect.txt \
  --k 5 \
  --out results/neighbors.tsv
```
### 4) Run everything end-to-end
```bash
libncker run \
  --gff /path/to/genome.gff \
  --lnc-ids-from-gff \
  --level transcript \
  --intersects /path/to/*_intersect.txt \
  --outdir results \
  --mode strict \
  --k 5
```
This now also writes:
- `intergenic.summary.tsv`
- `intergenic.distances.tsv`
- `intron.summary.tsv` if `--intron-tsv` is provided

### 5) Summarize intergenic distances from a GFF
```bash
libncker intergenic-stats \
  --gff /path/to/genome.gff.gz \
  --summary results/intergenic.summary.tsv \
  --distances-out results/intergenic.distances.tsv
```
### 6) Summarize intron lengths from a TSV column
```bash
libncker intron-stats \
  --tsv /path/to/introns.tsv \
  --column-index 8 \
  --summary results/intron.summary.tsv
```
### 7) Annotate DE mRNAs with NCBI gene summaries

Reads all `*_cis_regulation_module_output.txt` files from a directory, keeps only
differentially expressed (`DE_status == DE`) mRNAs, and fetches the NCBI gene summary
for each one. Produces one TSV per tissue.

```bash
libncker annotate \
  --cis-dir results/ \
  --outdir results/annotations/
```

With an NCBI API key (raises rate limit from 3 to 10 requests/sec):
```bash
libncker annotate \
  --cis-dir results/ \
  --outdir results/annotations/ \
  --ncbi-api-key YOUR_KEY
```

By default only first neighbors (`1_U` and `1_D`) are annotated, matching the intergenic
distance analysis. To include farther positions:
```bash
libncker annotate \
  --cis-dir results/ \
  --outdir results/annotations/ \
  --positions 1_U,1_D,2_U,2_D
```

#### eggNOG-mapper GO annotation (recommended for non-model organisms)

For organisms whose genes are not curated in NCBI (e.g. predicted *L. vannamei* RefSeq),
GO terms from NCBI Gene XML will be empty. Use `--emapper` to submit protein sequences to
the [eggNOG-mapper](http://eggnog-mapper.embl.de) web server and retrieve proper GO terms:

```bash
libncker annotate \
  --cis-dir results/ \
  --outdir results/annotations/ \
  --gff /path/to/genome.gff \
  --emapper
```

- `--gff` is optional but recommended: libncker reads the assembly accession from the GFF
  header, queries NCBI Taxonomy, and automatically selects the closest available eggNOG
  taxonomic scope (e.g. *Arthropoda* for shrimp).
- Without `--gff`, the scope defaults to *Arthropoda*.
- eggNOG-mapper GO terms overwrite the NCBI Gene XML GO columns when available.
- GO term names are resolved via the QuickGO API; raw GO IDs are used as fallback.

Output files: `<tissue>_de_mrna_annotations.tsv` — one file per tissue found, one row
per unique mRNA, with columns: `mRNA_ID`, `product`, `lncRNA_IDs`, `mRNA_positions`,
`ncbi_summary`, `go_biological_process`, `go_molecular_function`, `go_cellular_component`.

Works with any number of tissues; no configuration needed.

> **Note:** Uses only Python standard library — no extra packages required.
---

## Inputs
### Intersect tables (*_intersect.txt)

Each file should include:
- an ID column (lncRNA ID)
- an expression label column with values like:
  - `Up_<tissueA>_Down_<tissueB>`
  - `Down_<tissueA>_Up_<tissueB>`

Libncker infers tissue names from those labels and auto-detects the relevant columns.

### GFF (.gff/.gff3)

Used to:
- extract lncRNA IDs (`extract-ids` / `--lnc-ids-from-gff`)
- locate lncRNAs and mRNAs on contigs/chromosomes to compute neighbors
- compute intergenic distance statistics (`intergenic-stats`)
- optionally compute intron-length statistics during `run` with `--intron-tsv`

---

## Outputs
### In `--outdir`:
- `<tissue>_exclusive_lncRNAs_<mode>.txt`
- `<tissue>_up_consistent_but_incomplete.txt`
- `neighbors.tsv`
- `neighbors.summary.tsv`
- `intergenic.summary.tsv`
- `intergenic.distances.tsv`
- `intron.summary.tsv` if requested
- `summary.tsv`
- `<tissue>_cis_regulation_module_output.txt`

### From `annotate`:
- `<tissue>_de_mrna_annotations.tsv`

---

## How scoring works

Let N be the number of tissues detected from the intersect labels.

In strict mode, an lncRNA is “exclusive” to a tissue if it appears consistently as Up in that tissue across all N-1 relevant pairwise comparisons.

cis_score counts evidence across those comparisons; therefore the theoretical maximum in strict mode is:
- max(cis_score) = N - 1

For cross-dataset comparisons, consider reporting:
- cis_score_norm = cis_score / (N - 1)

---

## Common pitfalls
1) Do NOT include combined_intersect.txt together with pairwise intersects in strict mode.
Use: `*_intersect.txt`
2) Do not mix ID levels:
- `--level transcript` expects transcript IDs (e.g. XR_...)
- `--level gene` expects gene IDs (e.g. gene-LOC...)
3) If you lack some pairwise comparisons, strict mode will push many IDs into
`*_up_consistent_but_incomplete.txt`. Use `--mode lenient` `--min-required X`.
