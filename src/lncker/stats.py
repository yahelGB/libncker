from __future__ import annotations

import csv
import math
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, List

from lncker.gff_utils import iter_gff, norm_rna_id
from lncker.io_utils import iter_tsv_rows


def geometric_mean(values: Iterable[float]) -> float:
    cleaned = [float(v) for v in values if math.isfinite(float(v)) and float(v) > 0]
    if not cleaned:
        raise ValueError("No positive finite values available for geometric mean.")
    return math.exp(sum(math.log(v) for v in cleaned) / len(cleaned))


def median(values: List[float]) -> float:
    if not values:
        raise ValueError("Cannot compute median of an empty list.")
    xs = sorted(values)
    mid = len(xs) // 2
    if len(xs) % 2 == 1:
        return xs[mid]
    return (xs[mid - 1] + xs[mid]) / 2.0


@dataclass(frozen=True)
class GeneSpan:
    seqid: str
    start: int
    end: int
    strand: str
    gene_id: str


def _load_genes_from_gff(gff: Path) -> List[GeneSpan]:
    genes: List[GeneSpan] = []
    for feat in iter_gff(gff):
        if feat.ftype != "gene":
            continue
        gene_id = feat.attrs.get("ID", "")
        if not gene_id:
            continue
        genes.append(
            GeneSpan(
                seqid=feat.contig,
                start=feat.start,
                end=feat.end,
                strand=feat.strand,
                gene_id=norm_rna_id(gene_id),
            )
        )
    return genes


def _compute_intergenic_rows(genes: List[GeneSpan]) -> List[dict[str, object]]:
    grouped: dict[str, List[GeneSpan]] = {}
    for gene in genes:
        grouped.setdefault(gene.seqid, []).append(gene)

    out: List[dict[str, object]] = []
    for seqid, seq_genes in grouped.items():
        seq_genes.sort(key=lambda g: (g.start, g.end, g.gene_id))
        for i in range(len(seq_genes) - 1):
            a = seq_genes[i]
            b = seq_genes[i + 1]
            distance = b.start - a.end - 1
            if distance <= 0:
                continue
            out.append(
                {
                    "seqid": seqid,
                    "gene_A": a.gene_id,
                    "gene_B": b.gene_id,
                    "distance": distance,
                }
            )
    return out


def _write_tsv(path: Path, fieldnames: list[str], rows: Iterable[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def run_intergenic_stats(
    gff: Path,
    out_summary: Path | None = None,
    out_distances: Path | None = None,
) -> None:
    genes = _load_genes_from_gff(gff)
    if len(genes) < 2:
        raise ValueError("Not enough gene features were found in the GFF to compute intergenic distances.")

    rows = _compute_intergenic_rows(genes)
    if not rows:
        raise ValueError("No positive intergenic distances were computed.")

    distances = [float(row["distance"]) for row in rows]
    mean_bp = sum(distances) / len(distances)
    median_bp = median(distances)
    gmean_bp = geometric_mean(distances)

    print("=== Intergenic distances from GFF ===")
    print(f"N:\t{len(distances)}")
    print(f"Mean (bp):\t{mean_bp:.2f}")
    print(f"Median (bp):\t{median_bp:.2f}")
    print(f"Geometric mean (bp):\t{gmean_bp:.2f}")
    print(f"Mean (Kb):\t{mean_bp / 1000.0:.3f}")
    print(f"Median (Kb):\t{median_bp / 1000.0:.3f}")
    print(f"Geometric mean (Kb):\t{gmean_bp / 1000.0:.3f}")

    if out_summary is not None:
        _write_tsv(
            out_summary,
            ["metric", "value"],
            [
                {"metric": "input_gff", "value": str(gff)},
                {"metric": "n_intergenic_distances", "value": len(distances)},
                {"metric": "mean_bp", "value": f"{mean_bp:.6f}"},
                {"metric": "median_bp", "value": f"{median_bp:.6f}"},
                {"metric": "geometric_mean_bp", "value": f"{gmean_bp:.6f}"},
                {"metric": "mean_kb", "value": f"{mean_bp / 1000.0:.6f}"},
                {"metric": "median_kb", "value": f"{median_bp / 1000.0:.6f}"},
                {"metric": "geometric_mean_kb", "value": f"{gmean_bp / 1000.0:.6f}"},
            ],
        )

    if out_distances is not None:
        _write_tsv(out_distances, ["seqid", "gene_A", "gene_B", "distance"], rows)


def run_intron_stats(
    tsv: Path,
    column_index: int = 8,
    out_summary: Path | None = None,
) -> None:
    if column_index < 1:
        raise ValueError("column_index is 1-based and must be >= 1.")

    values: List[float] = []
    selected_idx = column_index - 1
    for row in iter_tsv_rows(tsv):
        if selected_idx >= len(row):
            continue
        raw = row[selected_idx].strip()
        try:
            value = float(raw)
        except ValueError:
            continue
        if math.isfinite(value) and value > 0:
            values.append(value)

    if not values:
        raise ValueError("No positive finite intron lengths found in the selected column.")

    median_nt = median(values)
    gmean_nt = geometric_mean(values)

    print("=== Intron lengths from TSV ===")
    print(f"N:\t{len(values)}")
    print(f"Median (nt):\t{median_nt:.2f}")
    print(f"Geometric mean (nt):\t{gmean_nt:.2f}")

    if out_summary is not None:
        _write_tsv(
            out_summary,
            ["metric", "value"],
            [
                {"metric": "input_tsv", "value": str(tsv)},
                {"metric": "column_index_1_based", "value": column_index},
                {"metric": "n_introns", "value": len(values)},
                {"metric": "median_nt", "value": f"{median_nt:.6f}"},
                {"metric": "geometric_mean_nt", "value": f"{gmean_nt:.6f}"},
            ],
        )
