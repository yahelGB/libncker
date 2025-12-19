from __future__ import annotations

import csv
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterator, List, Set, Tuple
from urllib.parse import unquote

from libncker.expr import parse_expression
from libncker.gff_utils import build_lnc_gene_to_transcripts, norm_rna_id, parse_attrs
from libncker.io_utils import iter_tsv_rows


@dataclass(frozen=True)
class Feature:
    contig: str
    ftype: str
    start: int
    end: int
    strand: str
    attrs: Dict[str, str]


def iter_gff_features(gff: Path) -> Iterator[Feature]:
    with gff.open("r", encoding="utf-8", errors="replace") as f:
        for line in f:
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9:
                continue
            contig, _src, ftype, start, end, _score, strand, _phase, attrs = parts
            yield Feature(
                contig=contig,
                ftype=ftype,
                start=int(start),
                end=int(end),
                strand=strand,
                attrs=parse_attrs(attrs),
            )


def load_exclusive_lncs(files: List[Path]) -> Dict[str, Set[str]]:
    out: Dict[str, Set[str]] = {}
    for p in files:
        tissue = p.name.split("_exclusive_lncRNAs", 1)[0]
        s: Set[str] = set()
        with p.open("r", encoding="utf-8", errors="replace") as f:
            for line in f:
                x = line.strip()
                if x:
                    s.add(x)
        out.setdefault(tissue, set()).update(s)
    return out


def _sort_key(f: Feature) -> Tuple[int, int, str]:
    return (f.start, f.end, f.ftype)


def distance_kb(lnc: Feature, mrna: Feature) -> float:
    if mrna.end < lnc.start:
        return (lnc.start - mrna.end) / 1000.0
    if mrna.start > lnc.end:
        return (mrna.start - lnc.end) / 1000.0
    overlap = min(lnc.end, mrna.end) - max(lnc.start, mrna.start) + 1
    return -overlap / 1000.0


def _build_lnc_index(contigs: Dict[str, List[Feature]]) -> Dict[str, List[Tuple[str, int]]]:
    idx: Dict[str, List[Tuple[str, int]]] = {}
    for contig_id, feats in contigs.items():
        for i, feat in enumerate(feats):
            if feat.ftype != "lnc_RNA":
                continue
            rid = norm_rna_id(feat.attrs.get("ID", ""))
            if rid:
                idx.setdefault(rid, []).append((contig_id, i))
    return idx


def find_neighbors_for_one_lnc(
    feats_sorted: List[Feature],
    lnc_index: int,
    k: int,
    exclusive_lncs_all: Set[str],
    target_lnc_id: str,
) -> Tuple[List[int], List[int]]:
    up_idxs: List[int] = []
    dn_idxs: List[int] = []
    used_parents_up: Set[str] = set()
    used_parents_dn: Set[str] = set()

    def is_blocking_exclusive_lnc(feat: Feature) -> bool:
        if feat.ftype != "lnc_RNA":
            return False
        rid = norm_rna_id(feat.attrs.get("ID", ""))
        return (rid in exclusive_lncs_all) and (rid != target_lnc_id)

    for i in range(lnc_index - 1, -1, -1):
        feat = feats_sorted[i]
        if is_blocking_exclusive_lnc(feat):
            break
        if feat.ftype != "mRNA":
            continue
        parent = feat.attrs.get("Parent", "")
        if not parent or parent in used_parents_up:
            continue
        used_parents_up.add(parent)
        up_idxs.append(i)
        if len(up_idxs) == k:
            break

    for i in range(lnc_index + 1, len(feats_sorted)):
        feat = feats_sorted[i]
        if is_blocking_exclusive_lnc(feat):
            break
        if feat.ftype != "mRNA":
            continue
        parent = feat.attrs.get("Parent", "")
        if not parent or parent in used_parents_dn:
            continue
        used_parents_dn.add(parent)
        dn_idxs.append(i)
        if len(dn_idxs) == k:
            break

    return up_idxs, dn_idxs


def _iter_intersect_rows(paths: List[Path]) -> Tuple[List[str], Iterator[List[str]]]:
    if not paths:
        raise ValueError("No intersect files provided.")
    header0 = next(iter_tsv_rows(paths[0]), None)
    if header0 is None:
        raise ValueError(f"First intersect file appears empty: {paths[0]}")

    def gen() -> Iterator[List[str]]:
        for p in paths:
            rows = iter_tsv_rows(p)
            _h = next(rows, None)
            for r in rows:
                if r:
                    yield r

    return header0, gen()


def _build_de_and_cis_counts_from_intersects(
    intersect_paths: List[Path],
) -> Tuple[Set[str], Set[str], Dict[str, Dict[str, int]], int, int]:
    header, rows = _iter_intersect_rows(intersect_paths)

    id_idx = header.index("ID") if "ID" in header else 0
    expr_idx = header.index("Expression") if "Expression" in header else (len(header) - 1)

    de_mrnas: Set[str] = set()
    tissues: Set[str] = set()

    seen_gene_comp: Set[Tuple[str, str, str]] = set()  # (gene, up, down) dedup
    mrna_up_counts: Dict[str, Dict[str, int]] = {}
    unique_comparisons: Set[Tuple[str, str]] = set()

    for r in rows:
        if id_idx >= len(r):
            continue
        gene_id = r[id_idx].strip()
        if not gene_id:
            continue

        expr = r[expr_idx].strip() if expr_idx < len(r) else ""
        parsed = parse_expression(expr)
        if parsed is None:
            continue

        up = parsed.up
        down = parsed.down

        de_mrnas.add(gene_id)
        tissues.add(up)
        tissues.add(down)
        unique_comparisons.add((up, down))

        key = (gene_id, up, down)
        if key in seen_gene_comp:
            continue
        seen_gene_comp.add(key)

        d = mrna_up_counts.setdefault(gene_id, {})
        d[up] = d.get(up, 0) + 1

    return de_mrnas, tissues, mrna_up_counts, len(tissues), len(unique_comparisons)


def run_neighbors(
    gff: Path,
    exclusive_files: List[Path],
    intersect_paths: List[Path],
    out_tsv: Path,
    k: int = 5,
    exclusive_level: str = "transcript",
) -> None:
    if exclusive_level not in {"transcript", "gene"}:
        raise ValueError("exclusive_level must be 'transcript' or 'gene'")

    tissue_to_ids = load_exclusive_lncs(exclusive_files)
    if not tissue_to_ids:
        raise ValueError("No exclusive lncRNA files provided or they are empty.")

    # If exclusive_level=gene, expand genes -> lnc_RNA transcripts using Parent links
    if exclusive_level == "gene":
        gene_to_tx = build_lnc_gene_to_transcripts(gff)

        tissue_to_lncs: Dict[str, Set[str]] = {}
        missing_genes = 0
        expanded_total = 0

        for tissue, genes in tissue_to_ids.items():
            txs: Set[str] = set()
            for gid in genes:
                # allow both "gene-LOC..." and "LOC..."
                keys = [gid]
                if not gid.startswith("gene-"):
                    keys.append("gene-" + gid)

                found_any = False
                for key in keys:
                    for t in gene_to_tx.get(key, set()):
                        txs.add(t)
                        found_any = True
                if not found_any:
                    missing_genes += 1
            tissue_to_lncs[tissue] = txs
            expanded_total += len(txs)

        if expanded_total == 0:
            raise ValueError(
                "Exclusive set is empty after level handling.\n"
                f"exclusive_level={exclusive_level}\n"
                "No lnc_RNA transcripts were found for your exclusive gene IDs in the GFF.\n"
                "Check that lnc_RNA lines have Parent=gene-... matching your gene IDs.\n"
                f"genes_without_any_lnc_RNA_children\t{missing_genes}\n"
            )
    else:
        tissue_to_lncs = tissue_to_ids

    exclusive_all: Set[str] = set()
    for s in tissue_to_lncs.values():
        exclusive_all |= s

    de_mrnas, tissues_inferred, mrna_up_counts, n_tissues, n_unique_comps = _build_de_and_cis_counts_from_intersects(
        intersect_paths
    )
    cis_expected = max(n_tissues - 1, 0)

    contigs: Dict[str, List[Feature]] = {}
    total_features = 0
    for feat in iter_gff_features(gff):
        contigs.setdefault(feat.contig, []).append(feat)
        total_features += 1

    for c in list(contigs.keys()):
        contigs[c].sort(key=_sort_key)

    lnc_index = _build_lnc_index(contigs)

    out_tsv.parent.mkdir(parents=True, exist_ok=True)
    summary_path = out_tsv.with_suffix(".summary.tsv")

    rows_written = 0
    lncs_total = sum(len(s) for s in tissue_to_lncs.values())
    lncs_found = 0
    lncs_missing = 0
    multi_hit_lncs = 0

    cis_rows_by_tissue: Dict[str, List[List[str]]] = {t: [] for t in tissue_to_lncs.keys()}

    with out_tsv.open("w", encoding="utf-8", newline="") as fo:
        w = csv.writer(fo, delimiter="\t")
        w.writerow([
            "tissue",
            "Contig",
            "lncRNA_ID",
            "lncRNA_strand",
            "mRNA_ID",
            "mRNA_strand",
            "mRNA_position",
            "distance_kb",
            "DE_status",
            "cis_score",
            "cis_expected",
            "cis_probability",
            "product",
        ])

        for tissue, lncs in tissue_to_lncs.items():
            for lnc_id in sorted(lncs):
                hits = lnc_index.get(lnc_id, [])
                if not hits:
                    lncs_missing += 1
                    continue
                lncs_found += 1
                if len(hits) > 1:
                    multi_hit_lncs += 1

                for contig_id, i in hits:
                    feats = contigs[contig_id]
                    lnc_feat = feats[i]

                    up_idxs, dn_idxs = find_neighbors_for_one_lnc(
                        feats_sorted=feats,
                        lnc_index=i,
                        k=k,
                        exclusive_lncs_all=exclusive_all,
                        target_lnc_id=lnc_id,
                    )

                    def emit(pos_label: str, mrna: Feature) -> None:
                        nonlocal rows_written
                        mrna_id = norm_rna_id(mrna.attrs.get("ID", "Unknown"))
                        product = unquote(mrna.attrs.get("product", "Unknown"))
                        de_status = "DE" if mrna_id in de_mrnas else "not_DE"

                        cis_score = mrna_up_counts.get(mrna_id, {}).get(tissue, 0)
                        cis_prob = (cis_score / cis_expected) if cis_expected > 0 else 0.0
                        dist = distance_kb(lnc_feat, mrna)

                        w.writerow([
                            tissue,
                            contig_id,
                            lnc_id,
                            lnc_feat.strand,
                            mrna_id,
                            mrna.strand,
                            pos_label,
                            f"{dist:.3f}",
                            de_status,
                            str(cis_score),
                            str(cis_expected),
                            f"{cis_prob:.3f}",
                            product,
                        ])
                        rows_written += 1

                        cis_rows_by_tissue[tissue].append([
                            contig_id,
                            lnc_id,
                            lnc_feat.strand,
                            mrna_id,
                            mrna.strand,
                            pos_label,
                            de_status,
                            str(cis_score),
                            f"{dist:.3f}",
                            product,
                        ])

                    for pos, idx in enumerate(up_idxs, start=1):
                        emit(f"{pos}_U", feats[idx])
                    for pos, idx in enumerate(dn_idxs, start=1):
                        emit(f"{pos}_D", feats[idx])

    for tissue, rows in cis_rows_by_tissue.items():
        out_cis = out_tsv.parent / f"{tissue}_cis_regulation_module_output.txt"
        with out_cis.open("w", encoding="utf-8", newline="") as fo:
            w = csv.writer(fo, delimiter="\t")
            w.writerow([
                "Contig",
                "lncRNA_ID",
                "lncRNA_strand",
                "mRNA_ID",
                "mRNA_strand",
                "mRNA_position",
                "DE_status",
                "cis_score",
                "distance_Kb",
                "product",
            ])
            for r in rows:
                w.writerow(r)

    with summary_path.open("w", encoding="utf-8") as s:
        s.write("key\tvalue\n")
        s.write(f"gff\t{gff}\n")
        s.write(f"exclusive_level\t{exclusive_level}\n")
        s.write(f"intersect_files_count\t{len(intersect_paths)}\n")
        s.write(f"unique_comparisons_detected\t{n_unique_comps}\n")
        s.write(f"exclusive_files_count\t{len(exclusive_files)}\n")
        s.write(f"exclusive_lncs_total\t{lncs_total}\n")
        s.write(f"exclusive_lncs_found_in_gff\t{lncs_found}\n")
        s.write(f"exclusive_lncs_missing_in_gff\t{lncs_missing}\n")
        s.write(f"exclusive_lncs_multi_hit_in_gff\t{multi_hit_lncs}\n")
        s.write(f"k_neighbors\t{k}\n")
        s.write(f"rows_written\t{rows_written}\n")
        s.write(f"contigs_loaded\t{len(contigs)}\n")
        s.write(f"features_loaded\t{total_features}\n")
        s.write(f"de_mrnas_in_intersects\t{len(de_mrnas)}\n")
        s.write(f"tissues_inferred_from_intersects\t{len(tissues_inferred)}\n")
        s.write(f"tissues_list\t{','.join(sorted(tissues_inferred))}\n")
        s.write(f"cis_expected\t{cis_expected}\n")