"""
deseq2.py — run all pairwise DESeq2 comparisons and produce lncker-compatible intersect tables.

Usage:
    python deseq2.py <count_table.txt> <metadata.csv> <output_dir> <padj> <log2fc> <min_counts>

Inputs:
    count_table.txt — genes x samples, tab-separated (output of merge_counts.py)
    metadata.csv    — columns: sample (row names), condition
    padj            — adjusted p-value cutoff (e.g. 0.05)
    log2fc          — minimum absolute log2 fold change (e.g. 1.5)
    min_counts      — minimum total counts across samples to retain a gene (e.g. 10)

Outputs per pairwise comparison (condA vs condB) written to output_dir:
    <condAvscondB>_intersect.txt        — ID + Expression columns (input for lncker)
    <condAvscondB>_results.txt          — full DESeq2 results sorted by padj
    <condAvscondB>_normalized_counts.txt — size-factor normalized counts for that pair
    <condAvscondB>_MA_plot.pdf          — MA plot

One shared output:
    PCA_plot.pdf                        — PCA on all samples (log1p raw counts)

Expression column format (IDEAMEX-compatible):
    Up_<condA>_Down_<condB>  when LFC > 0  (condA is higher)
    Down_<condA>_Up_<condB>  when LFC < 0  (condB is higher)
"""

import sys
import os
import itertools
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats

count_file  = sys.argv[1]
meta_file   = sys.argv[2]
output_dir  = sys.argv[3]
padj_cutoff = float(sys.argv[4])
lfc_cutoff  = float(sys.argv[5])
min_counts  = int(sys.argv[6])

os.makedirs(output_dir, exist_ok=True)

# --- load data ---
# count table is genes x samples → pydeseq2 expects samples x genes
counts   = pd.read_csv(count_file, sep="\t", index_col=0).T.astype(int)
metadata = pd.read_csv(meta_file, index_col=0)

common = counts.index.intersection(metadata.index)
if len(common) == 0:
    raise ValueError("No matching sample names between count table and metadata.")
counts   = counts.loc[common]
metadata = metadata.loc[common]

# --- PCA on all samples (log1p raw counts, no DESeq2 normalization needed) ---
log_counts = np.log1p(counts.values.astype(float))   # samples x genes
centered   = log_counts - log_counts.mean(axis=0)
_, S, Vt   = np.linalg.svd(centered, full_matrices=False)
coords     = centered @ Vt[:2].T
var_exp    = (S**2 / (S**2).sum())[:2] * 100

fig, ax = plt.subplots(figsize=(7, 5))
cond_col = metadata["condition"]
for cond, grp in cond_col.groupby(cond_col):
    idx = [counts.index.get_loc(s) for s in grp.index]
    ax.scatter(coords[idx, 0], coords[idx, 1], label=cond, s=80)
    for s, i in zip(grp.index, idx):
        ax.annotate(s, (coords[i, 0], coords[i, 1]), fontsize=7, ha="right")
ax.set_xlabel(f"PC1 ({var_exp[0]:.1f}%)")
ax.set_ylabel(f"PC2 ({var_exp[1]:.1f}%)")
ax.set_title("PCA — log1p raw counts (all samples)")
ax.legend()
plt.tight_layout()
plt.savefig(os.path.join(output_dir, "PCA_plot.pdf"))
plt.close()

# --- pairwise comparisons ---
conditions = sorted(metadata["condition"].unique())
pairs      = list(itertools.combinations(conditions, 2))

print(f"\nConditions detected: {conditions}")
print(f"Running {len(pairs)} pairwise comparison(s): {[f'{a}vs{b}' for a, b in pairs]}")
print(f"Thresholds: padj < {padj_cutoff}, |log2FC| >= {lfc_cutoff}, min_counts >= {min_counts}\n")

for cond_a, cond_b in pairs:
    tag = f"{cond_a}vs{cond_b}"
    print(f"{'='*50}")
    print(f"  {tag}")
    print(f"{'='*50}")

    # subset to just these two conditions
    mask       = metadata["condition"].isin([cond_a, cond_b])
    meta_sub   = metadata[mask].copy()
    counts_sub = counts.loc[meta_sub.index].copy()

    # pre-filter: keep genes with at least min_counts reads total in this subset
    counts_sub = counts_sub.loc[:, counts_sub.sum(axis=0) >= min_counts]

    # run DESeq2; contrast: cond_a vs cond_b → LFC = log2(cond_a / cond_b)
    dds = DeseqDataSet(counts=counts_sub, metadata=meta_sub, design="~condition")
    dds.deseq2()

    stat_res = DeseqStats(dds, contrast=["condition", cond_a, cond_b])
    stat_res.summary()
    results = stat_res.results_df.sort_values("padj", na_position="last")
    results.index.name = "ID"

    # full results
    results.to_csv(os.path.join(output_dir, f"{tag}_results.txt"), sep="\t")

    # normalized counts for this pair (genes x samples)
    norm = pd.DataFrame(
        dds.layers["normed_counts"],
        index=counts_sub.index,
        columns=counts_sub.columns,
    ).T
    norm.index.name = "ID"
    norm.to_csv(os.path.join(output_dir, f"{tag}_normalized_counts.txt"), sep="\t")

    # MA plot
    fig, ax = plt.subplots(figsize=(8, 5))
    sig  = results["padj"] < padj_cutoff
    base = np.log10(results["baseMean"] + 1)
    ax.scatter(base[~sig], results.loc[~sig, "log2FoldChange"],
               alpha=0.4, s=5, color="grey", label="Not significant")
    ax.scatter(base[sig],  results.loc[sig,  "log2FoldChange"],
               alpha=0.8, s=5, color="red",  label=f"padj < {padj_cutoff}")
    ax.axhline(0, color="blue", linewidth=0.8, linestyle="--")
    ax.set_xlabel("log10(mean expression + 1)")
    ax.set_ylabel("log2 fold change")
    ax.set_title(f"MA plot — {tag}")
    ax.legend(markerscale=3)
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, f"{tag}_MA_plot.pdf"))
    plt.close()

    # --- intersect table for lncker ---
    sig_mask = (
        results["padj"].notna() &
        (results["padj"] < padj_cutoff) &
        (results["log2FoldChange"].abs() >= lfc_cutoff)
    )
    sig_res = results[sig_mask].copy()

    def _expression_label(lfc):
        if lfc > 0:
            return f"Up_{cond_a}_Down_{cond_b}"
        return f"Down_{cond_a}_Up_{cond_b}"

    sig_res["Expression"] = sig_res["log2FoldChange"].apply(_expression_label)

    intersect = sig_res[["Expression"]].sort_index()
    intersect.index.name = "ID"
    intersect.to_csv(os.path.join(output_dir, f"{tag}_intersect.txt"), sep="\t")

    n_up = int((sig_res["log2FoldChange"] > 0).sum())
    n_dn = int((sig_res["log2FoldChange"] < 0).sum())
    print(f"  Significant (padj<{padj_cutoff}, |log2FC|>={lfc_cutoff}): {len(sig_res)}")
    print(f"    Up_{cond_a}_Down_{cond_b}: {n_up}")
    print(f"    Down_{cond_a}_Up_{cond_b}: {n_dn}")

print(f"\nDone. Output directory: {output_dir}")
