from __future__ import annotations

import argparse
from pathlib import Path

from libncker.annotate import run_annotate
from libncker.exclusive import run_exclusive
from libncker.gff_utils import extract_lncrna_ids_from_gff, write_extract_ids_summary
from libncker.neighbors import run_neighbors
from libncker.stats import run_intergenic_stats, run_intron_stats


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(prog="libncker", description="Libncker CLI")
    sub = p.add_subparsers(dest="cmd", required=True)

    # exclusive
    ex = sub.add_parser("exclusive", help="Compute tissue-exclusive lncRNAs from intersects tables.")
    ex.add_argument("--lnc-ids", required=True, help="Path to lncRNA IDs file (one ID per line).")
    ex.add_argument("--intersects", required=True, nargs="+", help="One or more *_intersect.txt files.")
    ex.add_argument("--outdir", required=True, help="Output directory.")
    ex.add_argument("--mode", choices=["strict", "lenient"], default="strict")
    ex.add_argument("--min-required", type=int, default=None)
    ex.add_argument("--level", choices=["transcript", "gene"], default="transcript")
    ex.add_argument("--gff", required=False, help="Required when --level gene (map transcript IDs -> gene parents).")

    # neighbors (cis module)
    nb = sub.add_parser("neighbors", help="Compute cis module outputs using exclusive lncRNAs + GFF + intersects.")
    nb.add_argument("--gff", required=True, help="Genome annotation GFF (NCBI/RefSeq style).")
    nb.add_argument("--exclusive", required=True, nargs="+", help="One or more *_exclusive_lncRNAs_*.txt files.")
    nb.add_argument("--intersects", required=True, nargs="+", help="One or more *_intersect.txt files.")
    nb.add_argument("--k", type=int, default=5, help="Neighbors K (default: 5).")
    nb.add_argument("--out", required=True, help="Output TSV path (e.g. results/neighbors.tsv).")
    nb.add_argument("--level", choices=["transcript", "gene"], default="transcript", help="Level of IDs in --exclusive files.")

    # extract-ids
    ei = sub.add_parser("extract-ids", help="Extract lncRNA IDs from GFF.")
    ei.add_argument("--gff", required=True, help="Genome annotation GFF.")
    ei.add_argument("--out", required=True, help="Output path (one ID per line).")
    ei.add_argument("--level", choices=["transcript", "gene"], default="transcript")
    ei.add_argument("--summary", required=False, help="Optional summary TSV path.")

    ig = sub.add_parser("intergenic-stats", help="Compute mean, median, and geometric mean intergenic distances from a GFF.")
    ig.add_argument("--gff", required=True, help="Genome annotation GFF.")
    ig.add_argument("--summary", required=False, help="Optional summary TSV path.")
    ig.add_argument("--distances-out", required=False, help="Optional detailed intergenic distances TSV path.")

    intron = sub.add_parser("intron-stats", help="Compute median and geometric mean intron lengths from a TSV column.")
    intron.add_argument("--tsv", required=True, help="Input TSV file.")
    intron.add_argument("--column-index", type=int, default=8, help="1-based column index with intron lengths (default: 8).")
    intron.add_argument("--summary", required=False, help="Optional summary TSV path.")

    # annotate
    an = sub.add_parser(
        "annotate",
        help="Annotate DE mRNAs from cis module outputs with NCBI gene summaries and Claude tissue-relevance annotations.",
    )
    an.add_argument("--cis-dir", required=True, help="Directory containing *_cis_regulation_module_output.txt files.")
    an.add_argument("--outdir", required=True, help="Output directory for per-tissue annotation TSV files.")
    an.add_argument("--ncbi-api-key", default="", help="NCBI API key (optional; raises rate limit from 3 to 10 req/s).")
    an.add_argument("--delay", type=float, default=0.4, help="Seconds between NCBI requests (default: 0.4).")

    # run = full pipeline
    runp = sub.add_parser("run", help="Run full pipeline: (IDs) -> exclusive -> cis module.")
    runp.add_argument("--gff", required=True, help="Genome annotation GFF.")
    group = runp.add_mutually_exclusive_group(required=True)
    group.add_argument("--lnc-ids", help="Path to lncRNA IDs file (one ID per line).")
    group.add_argument("--lnc-ids-from-gff", action="store_true", help="Extract lncRNA IDs from --gff.")
    runp.add_argument("--level", choices=["transcript", "gene"], default="transcript", help="Working level for exclusives.")
    runp.add_argument("--intersects", required=True, nargs="+", help="One or more *_intersect.txt files.")
    runp.add_argument("--outdir", required=True, help="Output directory for all results.")
    runp.add_argument("--mode", choices=["strict", "lenient"], default="strict")
    runp.add_argument("--min-required", type=int, default=None)
    runp.add_argument("--k", type=int, default=5, help="Neighbors K (default: 5).")
    runp.add_argument("--intron-tsv", required=False, help="Optional TSV for intron-length stats during the full run.")
    runp.add_argument(
        "--intron-column-index",
        type=int,
        default=8,
        help="1-based intron-length column for --intron-tsv (default: 8).",
    )

    return p


def main() -> None:
    parser = build_parser()
    args = parser.parse_args()

    if args.cmd == "annotate":
        run_annotate(
            cis_dir=Path(args.cis_dir),
            outdir=Path(args.outdir),
            ncbi_api_key=args.ncbi_api_key,
            delay=args.delay,
        )
        return

    if args.cmd == "exclusive":
        gff_path = Path(args.gff) if args.gff else None
        run_exclusive(
            intersect_paths=[Path(x) for x in args.intersects],
            lncrna_ids_path=Path(args.lnc_ids),
            outdir=Path(args.outdir),
            mode=args.mode,
            min_required=args.min_required,
            level=args.level,
            gff=gff_path,
        )
        return

    if args.cmd == "neighbors":
        run_neighbors(
            gff=Path(args.gff),
            exclusive_files=[Path(x) for x in args.exclusive],
            intersect_paths=[Path(x) for x in args.intersects],
            out_tsv=Path(args.out),
            k=args.k,
            exclusive_level=args.level,
        )
        return

    if args.cmd == "extract-ids":
        out = Path(args.out)
        out.parent.mkdir(parents=True, exist_ok=True)
        ids, stats = extract_lncrna_ids_from_gff(Path(args.gff), level=args.level)
        with out.open("w", encoding="utf-8") as f:
            for x in sorted(ids):
                f.write(x + "\n")
        if args.summary:
            write_extract_ids_summary(Path(args.summary), Path(args.gff), args.level, stats)
        return

    if args.cmd == "intergenic-stats":
        run_intergenic_stats(
            gff=Path(args.gff),
            out_summary=Path(args.summary) if args.summary else None,
            out_distances=Path(args.distances_out) if args.distances_out else None,
        )
        return

    if args.cmd == "intron-stats":
        run_intron_stats(
            tsv=Path(args.tsv),
            column_index=args.column_index,
            out_summary=Path(args.summary) if args.summary else None,
        )
        return

    if args.cmd == "run":
        outdir = Path(args.outdir)
        outdir.mkdir(parents=True, exist_ok=True)

        # 1) IDs
        if args.lnc_ids_from_gff:
            suffix = "tx" if args.level == "transcript" else "gene"
            lncs_path = outdir / f"lncRNA_IDs.from_gff.{suffix}.txt"
            ids, stats = extract_lncrna_ids_from_gff(Path(args.gff), level=args.level)
            with lncs_path.open("w", encoding="utf-8") as f:
                for x in sorted(ids):
                    f.write(x + "\n")
            write_extract_ids_summary(outdir / f"lncRNA_IDs.from_gff.{suffix}.summary.tsv", Path(args.gff), args.level, stats)
        else:
            lncs_path = Path(args.lnc_ids)

        # 2) exclusive
        run_exclusive(
            intersect_paths=[Path(x) for x in args.intersects],
            lncrna_ids_path=lncs_path,
            outdir=outdir,
            mode=args.mode,
            min_required=args.min_required,
            level=args.level,
            gff=Path(args.gff) if args.level == "gene" else None,
        )

        # 3) cis module (neighbors)
        exclusive_files = sorted(outdir.glob("*_exclusive_lncRNAs_*.txt"))
        neighbors_tsv = outdir / "neighbors.tsv"

        run_neighbors(
            gff=Path(args.gff),
            exclusive_files=exclusive_files,
            intersect_paths=[Path(x) for x in args.intersects],
            out_tsv=neighbors_tsv,
            k=args.k,
            exclusive_level=args.level,
        )

        run_intergenic_stats(
            gff=Path(args.gff),
            out_summary=outdir / "intergenic.summary.tsv",
            out_distances=outdir / "intergenic.distances.tsv",
        )

        if args.intron_tsv:
            run_intron_stats(
                tsv=Path(args.intron_tsv),
                column_index=args.intron_column_index,
                out_summary=outdir / "intron.summary.tsv",
            )
        return


if __name__ == "__main__":
    main()
