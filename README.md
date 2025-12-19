# libncker

Libncker is a lightweight CLI to:
1) identify tissue-exclusive lncRNAs from IDEAMEX intersect tables, and
2) compute a cis-neighborhood report (K nearest mRNAs) plus a simple `cis_score`.

---

## Table of contents
- [Install](#install)
- [Quickstart](#quickstart)
- [Commands](#commands)
- [Inputs](#inputs)
- [Outputs](#outputs)
- [How scoring works](#how-scoring-works)
- [Common pitfalls](#common-pitfalls)
- [Contributing](#contributing)
- [License](#license)

---

## Install

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
  --intersects /path/to/*vs*_intersect.txt \
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
  --intersects /path/to/*vs*_intersect.txt \
  --outdir results \
  --mode strict \
  --k 5
```

#### Note: if your GFF is compressed (.gff.gz), decompress first:
```bash
gunzip -c genome.gff.gz > genome.gff
```
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
### 2) Compute tissue-exclusive lncRNAs from IDEAMEX intersects
```bash
libncker exclusive \
  --lnc-ids lncRNA_IDs.txt \
  --level transcript \
  --intersects /path/to/*vs*_intersect.txt \
  --outdir results \
  --mode strict
```
Lenient mode (for incomplete pairwise sets or large N):
```bash
libncker exclusive \
  --lnc-ids lncRNA_IDs.txt \
  --level transcript \
  --intersects /path/to/*vs*_intersect.txt \
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
  --intersects /path/to/*vs*_intersect.txt \
  --outdir results \
  --mode strict \
  --k 5
```
---

## Inputs
### Intersect tables (*_intersect.txt)

Each file should include:
	•	an ID column (lncRNA ID)
	•	an expression label column with values like:
	•	Up_<tissueA>_Down_<tissueB>
	•	Down_<tissueA>_Up_<tissueB>

Libncker infers tissue names from those labels and auto-detects the relevant columns.

### GFF (.gff / .gff3)

Used to:
	•	extract lncRNA IDs (extract-ids / --lnc-ids-from-gff)
	•	locate lncRNAs and mRNAs on contigs to compute neighbors (neighbors)

---

## Outputs
### In --outdir:
- <tissue>_exclusive_lncRNAs_<mode>.txt
- <tissue>_up_consistent_but_incomplete.txt
- neighbors.tsv
- neighbors.summary.tsv
- summary.tsv
-<tissue>_cis_regulation_module_output.txt

---

## How scoring works

Let N be the number of tissues detected from the intersect labels.

In strict mode, an lncRNA is “exclusive” to a tissue if it appears consistently as Up in that tissue across all N-1 relevant pairwise comparisons.

cis_score counts evidence across those comparisons; therefore the theoretical maximum in strict mode is:
	•	max(cis_score) = N - 1

For cross-dataset comparisons, consider reporting:
	•	cis_score_norm = cis_score / (N - 1)

---

## Common pitfalls
	1.	Do NOT include combined_intersect.txt together with pairwise intersects in strict mode.
Use: *vs*_intersect.txt
	2.	Do not mix ID levels:
	•	--level transcript expects transcript IDs (e.g. XR_...)
	•	--level gene expects gene IDs (e.g. gene-LOC...)
	3.	If you lack some pairwise comparisons, strict mode will push many IDs into
*_up_consistent_but_incomplete.txt. Use --mode lenient --min-required X.
