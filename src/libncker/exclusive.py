from __future__ import annotations

from collections import Counter, defaultdict
from pathlib import Path
from typing import Optional

from libncker.expr import parse_expression
from libncker.gff_utils import build_transcript_to_gene_map
from libncker.io_utils import detect_id_and_expression_columns, iter_tsv_rows, read_ids_one_per_line


def _prefix(x: str) -> str:
    if x.startswith("gene-"):
        return "gene-"
    if "_" in x:
        return x.split("_", 1)[0]
    return x[:2]


def run_exclusive(
    intersect_paths: list[Path],
    lncrna_ids_path: Path,
    outdir: Path,
    mode: str = "strict",
    min_required: int | None = None,
    level: str = "transcript",
    gff: Optional[Path] = None,
) -> None:
    """
    Computes tissue-exclusive lncRNAs from IDEAMEX *_intersect.txt tables.

    level:
      - transcript: IDs are expected to match intersect ID column directly (XR_/NR_/etc)
      - gene: IDs list is gene-..., but intersects are transcript IDs -> we map transcript->Parent gene using --gff

    Output:
      - <tissue>_exclusive_lncRNAs_<mode>.txt (IDs at resolved level)
      - summary.tsv
    """
    if level not in {"transcript", "gene"}:
        raise ValueError("level must be 'transcript' or 'gene'")

    if not intersect_paths:
        raise ValueError("No intersect files provided.")

    lncs = read_ids_one_per_line(lncrna_ids_path)
    if not lncs:
        raise ValueError(f"lncRNA IDs file seems empty: {lncrna_ids_path}")

    if level == "gene" and gff is None:
        raise ValueError("When --level gene, you must provide --gff so we can map transcript IDs -> gene parents.")

    tx_to_gene = build_transcript_to_gene_map(gff) if (level == "gene" and gff is not None) else {}

    outdir.mkdir(parents=True, exist_ok=True)

    # per ID stats (ID here means resolved ID: transcript or gene)
    total = defaultdict(int)     # id -> total occurrences (filtered)
    up_counts = defaultdict(int) # (id, up_tissue) -> count
    up_set = defaultdict(set)    # id -> set(up_tissue)

    tissues = set()
    comparisons = set()

    invalid_expr_count = 0
    kept_rows = 0
    skipped_not_lnc = 0
    skipped_unmapped_to_gene = 0
    files_processed = 0

    # debug prefix tallies
    prefixes_ids = Counter(_prefix(x) for x in lncs)
    prefixes_intersects = Counter()

    for p in intersect_paths:
        rows = iter_tsv_rows(p)
        try:
            header = next(rows)
        except StopIteration:
            continue

        cols = detect_id_and_expression_columns(header)
        if cols is None:
            raise ValueError(f"Intersect file missing ID/Expression header: {p}")

        id_idx, expr_idx = cols
        files_processed += 1

        for r in rows:
            if not r:
                continue
            if id_idx >= len(r) or expr_idx >= len(r):
                continue

            raw_id = r[id_idx].strip()
            expr = r[expr_idx].strip()

            if not raw_id or raw_id == "ID":
                continue

            prefixes_intersects[_prefix(raw_id)] += 1

            parsed = parse_expression(expr)
            if parsed is None:
                invalid_expr_count += 1
                continue

            tissues.add(parsed.up)
            tissues.add(parsed.down)
            comparisons.add(frozenset((parsed.up, parsed.down)))

            if level == "transcript":
                resolved_id = raw_id
            else:
                # level == gene: map transcript ID -> Parent gene-...
                resolved_id = tx_to_gene.get(raw_id, "")
                if not resolved_id:
                    skipped_unmapped_to_gene += 1
                    continue

            if resolved_id not in lncs:
                skipped_not_lnc += 1
                continue

            kept_rows += 1
            total[resolved_id] += 1
            up_counts[(resolved_id, parsed.up)] += 1
            up_set[resolved_id].add(parsed.up)

    if kept_rows == 0:
        # same style of helpful error you liked earlier
        raise ValueError(
            "No rows from intersect files matched your --lnc-ids list.\n"
            f"lnc-ids file: {lncrna_ids_path}\n"
            "Example: if your intersects contain XM_/NM_ (mRNAs) but your lnc list is XR_/NR_ (lncRNAs), exclusives will be 0.\n"
            f"level={level}\n"
            f"Top prefixes in lnc-ids: {dict(prefixes_ids.most_common(10))}\n"
            f"Top prefixes in intersects: {dict(prefixes_intersects.most_common(10))}\n"
            + (f"unmapped_transcripts_to_gene: {skipped_unmapped_to_gene}\n" if level == "gene" else "")
        )

    N = len(tissues)
    expected_comparisons = (N * (N - 1)) // 2
    observed_comparisons = len(comparisons)

    required_strict = N - 1
    if mode == "strict":
        required = required_strict
    else:
        required = min_required if (min_required is not None) else required_strict

    exclusives: dict[str, list[str]] = {t: [] for t in sorted(tissues)}
    mixed_ids = 0
    up_consistent_incomplete: dict[str, list[str]] = {t: [] for t in sorted(tissues)}  # strict only

    for resolved_id, cnt in total.items():
        ups = up_set[resolved_id]
        if len(ups) != 1:
            mixed_ids += 1
            continue
        up_t = next(iter(ups))

        if mode == "strict":
            if cnt == required and up_counts[(resolved_id, up_t)] == required:
                exclusives[up_t].append(resolved_id)
            else:
                up_consistent_incomplete[up_t].append(resolved_id)
        else:
            if cnt >= required and up_counts[(resolved_id, up_t)] == cnt:
                exclusives[up_t].append(resolved_id)

    # write outputs
    for t, ids in exclusives.items():
        ids_sorted = sorted(set(ids))
        out = outdir / f"{t}_exclusive_lncRNAs_{mode}.txt"
        with out.open("w", encoding="utf-8") as f:
            for x in ids_sorted:
                f.write(x + "\n")

    if mode == "strict":
        for t, ids in up_consistent_incomplete.items():
            ids_sorted = sorted(set(ids))
            out = outdir / f"{t}_up_consistent_but_incomplete.txt"
            with out.open("w", encoding="utf-8") as f:
                for x in ids_sorted:
                    f.write(x + "\n")

    # summary
    summary_path = outdir / "summary.tsv"
    with summary_path.open("w", encoding="utf-8") as s:
        s.write("key\tvalue\n")
        s.write(f"exclusive_level\t{level}\n")
        s.write(f"lncrna_ids_count\t{len(lncs)}\n")
        s.write(f"files_processed\t{files_processed}\n")
        s.write(f"kept_rows_lncrna\t{kept_rows}\n")
        s.write(f"skipped_not_lnc_rows\t{skipped_not_lnc}\n")
        s.write(f"invalid_expression_rows\t{invalid_expr_count}\n")
        if level == "gene":
            s.write(f"unmapped_transcripts_to_gene\t{skipped_unmapped_to_gene}\n")
        s.write(f"tissues_inferred\t{N}\n")
        s.write(f"tissues_list\t{','.join(sorted(tissues))}\n")
        s.write(f"expected_pairwise_comparisons\t{expected_comparisons}\n")
        s.write(f"observed_pairwise_comparisons\t{observed_comparisons}\n")
        s.write(f"required_patterns_strict\t{required_strict}\n")
        s.write(f"mode\t{mode}\n")
        s.write(f"required_patterns_used\t{required}\n")
        s.write(f"ids_mixed_up_direction\t{mixed_ids}\n")
        for t in sorted(tissues):
            out = outdir / f"{t}_exclusive_lncRNAs_{mode}.txt"
            try:
                n = sum(1 for _ in out.open("r", encoding="utf-8"))
            except FileNotFoundError:
                n = 0
            s.write(f"exclusive_count_{t}\t{n}\n")