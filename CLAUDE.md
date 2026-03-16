# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

**libncker** (Lightweight lncRNA & Cis-neighborhood Analysis) is a pure Python 3.9+ CLI bioinformatics tool with zero external dependencies. It analyzes long non-coding RNAs (lncRNAs) from RNA-seq differential expression data through three core functions:

1. **Tissue-exclusive lncRNA identification** — finds lncRNAs uniquely upregulated in specific tissues
2. **Cis-neighborhood scoring** — identifies K nearest mRNA neighbors and computes cis-regulatory scores
3. **Intergenic/intron distance statistics** — characterizes genomic distances using geometric means and medians

## Commands

```bash
# Install for development
pip install -e .

# Run CLI
libncker -h
python3 -m libncker -h

# Run full pipeline
libncker run --gff <genome.gff.gz> --intersects <dir/> --output <dir/>

# Individual subcommands
libncker extract-ids --gff <genome.gff.gz> --output <dir/>
libncker exclusive --ids <lnc_ids.txt> --intersects <dir/> --output <dir/>
libncker neighbors --gff <genome.gff.gz> --exclusive-dir <dir/> --intersects <dir/> --output <dir/>
libncker intergenic-stats --gff <genome.gff.gz> --output <dir/>
libncker intron-stats --tsv <file.tsv> --column 3 --output <dir/>

# Build distribution
python3 -m build
```

No test suite or linter is currently configured.

## Architecture

### Pipeline Data Flow

```
GFF (genome annotation)
    │
    ├──► extract-ids  ──► lnc_ids.txt
    │
    ├──► exclusive (+ lnc_ids.txt + intersect tables) ──► <tissue>_exclusive_lncRNAs_<mode>.txt
    │
    ├──► neighbors (+ exclusive files + intersect tables) ──► neighbors.tsv, cis_regulation_module_output.txt
    │
    └──► intergenic-stats ──► intergenic.summary.tsv
```

### Module Responsibilities

| Module | Role |
|--------|------|
| `__main__.py` | CLI parser (`argparse`) and command dispatcher |
| `gff_utils.py` | GFF3 parsing; `GffFeature` frozen dataclass; ID/parent extraction; supports `.gz` |
| `io_utils.py` | File I/O: `open_text()` for gz/plain, TSV iteration, column detection |
| `expr.py` | Parses IDEAMEX-style expression labels (`Up_<tissue>_Down_<other>`) into `ParsedExpression` |
| `exclusive.py` | Core logic to identify tissue-exclusive lncRNAs (strict and lenient modes) |
| `neighbors.py` | K-nearest mRNA neighbor search on sorted GFF index; computes cis_probability scores |
| `stats.py` | Geometric mean and median computations for intergenic distances and intron lengths |

### Key Design Details

**ID levels:** Two modes — `transcript` (XR_, NM_) and `gene` (gene-LOC...). The `neighbors` module can expand gene-level IDs to transcripts via GFF `Parent` links at runtime.

**Exclusive modes:**
- **Strict:** lncRNA appears "up" in same tissue across all (N−1) pairwise comparisons
- **Lenient:** lncRNA appears "up" ≥ `min_required` times, all in the same tissue

**GFF indexing in `neighbors.py`:** Features are grouped by contig and sorted by start position. Neighbor search is bidirectional (upstream/downstream), stopping at other exclusive lncRNAs.

**Cis probability:** `cis_score / cis_expected`, where `cis_expected = max_tissues − 1` (from the intersect tables passed to neighbors).

**No external dependencies:** Only Python stdlib (`argparse`, `pathlib`, `dataclasses`, `csv`, `gzip`, `math`).
