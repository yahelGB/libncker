from __future__ import annotations

import argparse
from pathlib import Path

from libncker.exclusive import run_exclusive
from libncker.gff_utils import extract_lncrna_ids_from_gff, write_extract_ids_summary
from libncker.neighbors import run_neighbors


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(prog="libncker", description="Libncker CLI")
    sub = p.add_subparsers(dest="cmd", required=True)

    # exclusive
    ex = sub.add_parser("exclusive", help="Compute tissue-exclusive lncRNAs from IDEAMEX intersects.")
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

    return p


def main() -> None:
    parser = build_parser()
    args = parser.parse_args()

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
        return


if __name__ == "__main__":
    main()