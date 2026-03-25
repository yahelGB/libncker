"""
merge_counts.py — merge per-sample coverage files into a genes x samples count table.

Usage:
    python merge_counts.py <output_file> <sample1.coverage> [<sample2.coverage> ...]

Input format (per sample, from: samtools coverage | awk '{print $1"\t"$4}'):
    #rname    numreads
    gene1     150
    gene2     0
    ...

Output: tab-separated count matrix (genes x samples), zero-sum rows removed.
"""

import sys
import os
import pandas as pd

output_file = sys.argv[1]
input_files = sys.argv[2:]

dfs = []
for f in input_files:
    sample = os.path.basename(f).replace(".coverage", "")
    # comment='#' skips the "#rname numreads" header line
    df = pd.read_csv(
        f,
        sep="\t",
        comment="#",
        header=None,
        names=["gene", "numreads"]
    )
    df = df.set_index("gene")[["numreads"]].rename(columns={"numreads": sample})
    dfs.append(df)

table = pd.concat(dfs, axis=1).fillna(0).astype(int)
table = table[table.sum(axis=1) > 0]

os.makedirs(os.path.dirname(os.path.abspath(output_file)), exist_ok=True)
table.to_csv(output_file, sep="\t")

print(f"Count table written: {table.shape[0]} genes x {table.shape[1]} samples → {output_file}")
