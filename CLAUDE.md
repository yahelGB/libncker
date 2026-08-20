# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

**lncker** (Lightweight lncRNA & Cis-neighborhood Analysis) is a pure Python 3.9+ CLI bioinformatics tool with **zero required dependencies** (stdlib only; the `annotate` command reaches out to NCBI/eggNOG-mapper/QuickGO web APIs at runtime but still uses only `urllib`). It analyzes long non-coding RNAs (lncRNAs) from RNA-seq differential expression data through four core functions:

1. **Tissue-exclusive lncRNA identification** — finds lncRNAs uniquely upregulated in specific tissues
2. **Cis-neighborhood scoring** — identifies K nearest mRNA neighbors and computes cis-regulatory scores
3. **Intergenic/intron distance statistics** — characterizes genomic distances using geometric means and medians
4. **DE mRNA annotation** — annotates cis-neighbor mRNAs with NCBI gene summaries and GO terms

An upstream Snakemake RNA-seq pipeline (`snakemake/`) produces the `*_intersect.txt` inputs; it has its own `snakemake/CLAUDE.md`.

## Commands

```bash
# Install for development
pip install -e .

# Run CLI
lncker -h
python3 -m lncker -h

# Run without installing (no deps, so this always works)
PYTHONPATH=src python3 -m lncker -h

# Run full pipeline (extract-ids -> exclusive -> neighbors -> intergenic-stats)
lncker run --gff <genome.gff[.gz]> --lnc-ids-from-gff --level transcript \
  --intersects <dir>/*_intersect.txt --outdir <dir> --mode strict --k 5

# Individual subcommands (note: --intersects takes globbed FILES, not a directory)
lncker extract-ids     --gff <genome.gff> --level transcript --out lnc_ids.txt [--summary s.tsv]
lncker exclusive       --lnc-ids lnc_ids.txt --intersects <dir>/*_intersect.txt --outdir <dir> --mode strict
lncker neighbors       --gff <genome.gff> --exclusive <dir>/*_exclusive_lncRNAs_*.txt --intersects <dir>/*_intersect.txt --k 5 --out neighbors.tsv
lncker intergenic-stats --gff <genome.gff> [--summary s.tsv] [--distances-out d.tsv]
lncker intron-stats    --tsv <file.tsv> --column-index 8 [--summary s.tsv]
lncker annotate        --cis-dir <dir> --outdir <dir> [--emapper --gff <genome.gff>]

# Build distribution
python3 -m build
```

### Regression check (there is no test suite or linter)

`examples/output_test/` is the de-facto golden fixture, and it was generated at **`--level gene`** — not transcript. Reproduce it with:

```bash
PYTHONPATH=src python3 -m lncker run \
  --gff GCF_003789085.1_ASM378908v1_genomic.gff \
  --lnc-ids-from-gff --level gene \
  --intersects examples/data_test/*_intersect.txt \
  --outdir /tmp/regress --mode strict --k 5
diff -r /tmp/regress examples/output_test
```

Takes ~10 s. Every file should match byte-for-byte **except**:
- the three `*summary*.tsv` files, which embed the absolute `--gff` path (author's macOS path in the fixture) — diff only line 2;
- `lncRNA_IDs.from_gff.tx.txt` / `.tx.summary.tsv` in the fixture, which use the pre-`cc8c57c` filename suffix; current code writes `.gene.` / `.transcript.`.

Caveats when working from a fresh clone:
- `.gitignore` excludes `*.gff`, so the ~191 MB `GCF_003789085.1_ASM378908v1_genomic.gff` (*Litopenaeus vannamei* RefSeq) at the repo root is **untracked and local-only**, and `examples/data_test/` ships no GFF. Without it only `exclusive` and `intron-stats` can run; `neighbors`, `intergenic-stats`, and `run` cannot.
- Running `exclusive` against `examples/data_test/lncRNA_IDs.txt` (transcript-level) will **not** reproduce `output_test` — that path yields transcript IDs and slightly different per-tissue counts. Use the gene-level `run` above.

`examples/*.R` are downstream ggplot2 analysis scripts (cis_score vs. distance, first-neighbor distances, annotation plots) that consume `output_test/`. They are not part of the Python package and have hardcoded author paths.

## Architecture

### Pipeline Data Flow

```
GFF (genome annotation)
    │
    ├──► extract-ids  ──► lncRNA_IDs.from_gff.<level>.txt
    │
    ├──► exclusive (+ lnc IDs + intersect tables) ──► <tissue>_exclusive_lncRNAs_<mode>.txt
    │                                                  <tissue>_up_consistent_but_incomplete.txt
    │
    ├──► neighbors (+ exclusive files + intersect tables) ──► neighbors.tsv, neighbors.summary.tsv,
    │                                                          <tissue>_cis_regulation_module_output.txt
    │
    ├──► intergenic-stats ──► intergenic.summary.tsv, intergenic.distances.tsv
    │
    └──► annotate (reads *_cis_regulation_module_output.txt) ──► <tissue>_de_mrna_annotations.tsv
```

`run` chains extract-ids → exclusive → neighbors → intergenic-stats (and intron-stats if `--intron-tsv` given). `annotate` is a separate, network-dependent step run after `run`.

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
| `annotate.py` | Reads cis module outputs, fetches NCBI gene summaries/GO terms for DE mRNAs (`urllib` → NCBI eutils, QuickGO) |
| `emapper.py` | eggNOG-mapper GO annotation: web API (`eggnog-mapper.embl.de`) or local `emapper.py`; taxon auto-detection from GFF via NCBI Taxonomy |

### Key Design Details

**ID levels:** Two modes — `transcript` (XR_, NM_) and `gene` (gene-LOC...). The `neighbors` and `exclusive` modules can expand gene-level IDs to transcripts via GFF `Parent` links at runtime (`--gff` required when `--level gene` for `exclusive`).

**Exclusive modes:**
- **Strict:** lncRNA appears "up" in same tissue across all (N−1) pairwise comparisons
- **Lenient:** lncRNA appears "up" ≥ `--min-required` times, all in the same tissue

**Hard dependency on NCBI/RefSeq GFF conventions.** All GFF logic keys off exact type and attribute strings, so a GFF from another source generally will not work without adaptation:
- feature types must be literally `gene`, `lnc_RNA`, and `mRNA` (Ensembl-style `transcript`/`ncRNA` yields zero hits, and the failure is silent — empty outputs, not an error);
- lncRNA discovery uses `gene_biotype=lncRNA` on genes **unioned with** the `Parent` of every `lnc_RNA` feature, so it still works when `gene_biotype` is absent;
- transcript IDs are normalized by stripping a leading `rna-` (`norm_rna_id`), gene IDs keep their `gene-` prefix.

**GFF indexing in `neighbors.py`:** Features are grouped by contig and sorted by `(start, end, ftype)`. Neighbor search walks outward from the lncRNA's index in both directions, collecting `mRNA` features and **deduplicating by `Parent`** so each neighboring gene contributes one transcript. The walk stops early at any *other* exclusive lncRNA (the target itself does not block). mRNA positions are labelled `<k>_U` / `<k>_D` (upstream/downstream); `annotate` defaults to `1_U,1_D` (first neighbors). `distance_kb` is **negative when lncRNA and mRNA overlap** (it reports `-overlap/1000`), which downstream consumers must handle.

**Cis score / probability:** `cis_score` counts supporting pairwise comparisons; in strict mode `max(cis_score) = N − 1` where N is the number of tissues detected from intersect labels. `neighbors.tsv` already emits `cis_expected` (= N − 1) and `cis_probability` (= `cis_score / cis_expected`); the per-tissue `*_cis_regulation_module_output.txt` drops those two columns and carries the raw `cis_score` only.

**Intersect header contract:** files must have literal `ID` and `Expression` column headers, but the two consumers disagree on strictness — `io_utils.detect_id_and_expression_columns` (used by `exclusive`) raises a clear error when either is missing, while `neighbors._build_de_and_cis_counts_from_intersects` silently falls back to column 0 and the last column. A malformed header therefore fails loudly in `exclusive` but produces quiet garbage in `neighbors`. `neighbors` also takes the header from the **first** file only and assumes all files share that layout.

**Diagnostic summaries over exceptions:** every command writes a `key\tvalue` summary TSV (`summary.tsv`, `neighbors.summary.tsv`, `intergenic.summary.tsv`, `lncRNA_IDs.from_gff.<level>.summary.tsv`) recording rows kept/skipped, tissues inferred, unmapped IDs, and expected-vs-observed comparison counts. When output looks wrong, read those counters first — most silent mismatches (ID-level confusion, missing comparisons) show up as `skipped_not_lnc_rows`, `unmapped_transcripts_to_gene`, or `exclusive_lncs_missing_in_gff`. The `exclusive` and `neighbors` modules do raise `ValueError` with prefix tallies when *nothing* matches at all.

**`annotate` networking:** Uses NCBI eutils (rate-limited: 3 req/s, or 10 with `--ncbi-api-key`, throttled by `--delay`). For non-model organisms whose NCBI GO terms are empty, `--emapper` submits proteins to eggNOG-mapper (web) or runs a local DB (`--emapper-db`); taxonomic scope auto-detected from `--gff` header or set via `--tax-scope`. Still stdlib-only — no pip packages.

**No external Python dependencies:** stdlib only (`argparse`, `pathlib`, `dataclasses`, `csv`, `gzip`, `math`, `urllib`, `json`).

## Common Pitfalls

- `--intersects` expects **globbed pairwise files** (`*_intersect.txt`), not a directory and not a combined intersect file mixed with pairwise ones.
- Do not mix ID levels: `--level transcript` expects `XR_...`; `--level gene` expects `gene-LOC...`.
- Missing pairwise comparisons push IDs into `*_up_consistent_but_incomplete.txt` under strict mode — use `--mode lenient --min-required X`. Those `*_up_consistent_but_incomplete.txt` files are written **only in strict mode**; lenient mode silently drops the near-misses.
- `neighbors --exclusive` globs `*_exclusive_lncRNAs_*.txt` and derives the tissue name from the filename prefix (text before `_exclusive_lncRNAs`) — renaming those files reassigns tissues.
- `--level` must be passed consistently to `exclusive` and `neighbors`; `run` handles this for you and only forwards `--gff` to `exclusive` when `--level gene`.
