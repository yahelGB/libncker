# `examples/output_test/` — regenerating this directory

Every file here (except `annotations/`) is produced by one command and is
byte-reproducible from the committed inputs **plus** the pinned reference GFF.

## The pinned input

The reference annotation is **not committed** — it is ~191 MB and `.gitignore`
excludes `*.gff`. Download it into the repository root:

```bash
curl -O https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/003/789/085/GCF_003789085.1_ASM378908v1/GCF_003789085.1_ASM378908v1_genomic.gff.gz
gunzip GCF_003789085.1_ASM378908v1_genomic.gff.gz
```

| | |
|---|---|
| accession | `NCBI_Assembly:GCF_003789085.1` |
| assembly | `ASM378908v1` (*Penaeus vannamei*) |
| annotation release | NCBI *Penaeus vannamei* Annotation Release 100 |
| md5 | `2b9e422c5ef67de99eaefe6b74d09b31` |
| bytes | `191504969` |

`PROVENANCE.txt` in the repository root records the same values and is written by
`lncker run --write-gff-provenance`.

## The command

Run from the repository root. Note `--level gene` — the outputs here are
gene-level, and `--level transcript` produces different per-tissue counts.

```bash
lncker run \
  --gff GCF_003789085.1_ASM378908v1_genomic.gff \
  --lnc-ids-from-gff --level gene \
  --intersects examples/data_test/*_intersect.txt \
  --outdir examples/output_test \
  --mode strict --k 5
```

Takes about 10 s. Everything except `annotations/` is overwritten.

## Verifying a regeneration

Each file opens with a `#` provenance header recording the tool version, git
SHA, input md5s, and the full command line, so the header differs on every run
by design. Compare the data only:

```bash
git stash && lncker run ... --outdir /tmp/check   # command above
for f in examples/output_test/*.tsv examples/output_test/*.txt; do
  diff <(grep -v '^#' "$f") <(grep -v '^#' "/tmp/check/$(basename "$f")") \
    && echo "OK $(basename "$f")"
done
```

Read these files with the comment character set:

```r
read.delim("neighbors.tsv", comment.char = "#")          # R
```
```python
pd.read_csv("neighbors.tsv", sep="\t", comment="#")      # pandas
```

## `annotations/` is different

`annotations/*_de_mrna_annotations.tsv` come from `lncker annotate`, which
queries NCBI eutils and QuickGO live. They depend on the state of those
databases on the day they were fetched and are **not** reproducible offline or
byte-for-byte. They are kept as illustrative output only. To refresh them:

```bash
lncker annotate --cis-dir examples/output_test --outdir examples/output_test/annotations
```
