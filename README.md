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
