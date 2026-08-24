from __future__ import annotations

import argparse
from pathlib import Path

import sys

from lncker.annotate import run_annotate
from lncker.exclusive import run_exclusive
from lncker.gff_utils import extract_lncrna_ids_from_gff, write_extract_ids_summary
from lncker.neighbors import run_neighbors
from lncker.provenance import Provenance, write_gff_provenance_file
from lncker.stats import run_intergenic_stats, run_intron_stats
from lncker.validate import (
    DEFAULT_MIN_ID_OVERLAP,
    DropLog,
    InputContractError,
    check_expression_ids_in_gff,
    set_default_drop_log,
    validate_gff,
)


def _contract_flags(sp: argparse.ArgumentParser) -> None:
    """Flags controlling input-contract enforcement (see lncker.validate)."""
    sp.add_argument(
        "--no-validate",
        action="store_true",
        default=False,
        help="Skip input-contract validation. Use only when you know the annotation is non-RefSeq.",
    )
    sp.add_argument(
        "--min-id-overlap",
        type=float,
        default=DEFAULT_MIN_ID_OVERLAP,
        help=(
            "Minimum fraction of expression-matrix IDs that must be present in the GFF "
            f"(default: {DEFAULT_MIN_ID_OVERLAP}). Below this, lncker aborts."
        ),
    )


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(prog="lncker", description="Lncker CLI")
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
        help="Annotate DE mRNAs from cis module outputs with NCBI gene summaries and GO terms.",
    )
    an.add_argument("--cis-dir", required=True, help="Directory containing *_cis_regulation_module_output.txt files.")
    an.add_argument("--outdir", required=True, help="Output directory for per-tissue annotation TSV files.")
    an.add_argument("--positions", default="1_U,1_D", help="Comma-separated mRNA positions to include (default: 1_U,1_D).")
    an.add_argument("--ncbi-api-key", default="", help="NCBI API key (optional; raises rate limit from 3 to 10 req/s).")
    an.add_argument("--delay", type=float, default=0.4, help="Seconds between NCBI requests (default: 0.4).")
    an.add_argument(
        "--emapper",
        action="store_true",
        default=False,
        help=(
            "Run eggNOG-mapper for GO annotation (recommended for non-model organisms). "
            "Requires internet access to eggnog-mapper.embl.de. "
            "Auto-detects the organism's taxon from --gff if provided."
        ),
    )
    an.add_argument(
        "--gff",
        default=None,
        help="Genome annotation GFF (used with --emapper to auto-detect the organism's taxonomic scope).",
    )
    an.add_argument(
        "--emapper-db",
        default=None,
        help=(
            "Path to local eggNOG database directory. "
            "Triggers local mode: runs emapper.py on PATH instead of the web API. "
            "Useful on HPC nodes with no internet or when the DB is already downloaded."
        ),
    )
    an.add_argument(
        "--tax-scope",
        type=int,
        default=None,
        help=(
            "Explicit eggNOG taxonomic scope ID (e.g. 6656=Arthropoda, 33208=Metazoa). "
            "Overrides auto-detection from --gff. Useful when NCBI is unreachable."
        ),
    )
    an.add_argument(
        "--emapper-cpu",
        type=int,
        default=4,
        help="CPUs for local emapper.py run (default: 4, ignored in web mode).",
    )

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
    runp.add_argument(
        "--write-gff-provenance",
        action="store_true",
        default=False,
        help="Also write PROVENANCE.txt next to --gff (accession, release, md5).",
    )

    for sp in (ex, nb, ei, ig, intron, an, runp):
        _contract_flags(sp)

    return p


def _intersect_ids(paths: list[Path]) -> set[str]:
    """Union of the ID column across intersect files, for the contract check."""
    from lncker.io_utils import detect_id_and_expression_columns, iter_tsv_rows

    ids: set[str] = set()
    for path in paths:
        rows = iter_tsv_rows(path)
        try:
            header = next(rows)
        except StopIteration:
            raise InputContractError("intersect_header", f"{path} is empty")
        cols = detect_id_and_expression_columns(header)
        if cols is None:
            raise InputContractError(
                "intersect_header",
                f"{path} must have literal 'ID' and 'Expression' column headers; got {header[:8]}",
            )
        id_idx, _ = cols
        for row in rows:
            if id_idx < len(row):
                value = row[id_idx].strip()
                if value:
                    ids.add(value)
    return ids


def _enforce_contract(args: argparse.Namespace, drops: DropLog) -> None:
    """Validate the declared input contract before any work is done.

    Skipped entirely under --no-validate. Any failure raises InputContractError
    naming the check that failed, which main() turns into exit code 2.
    """
    if getattr(args, "no_validate", False):
        print("[contract] validation skipped (--no-validate)", file=sys.stderr)
        return

    gff = getattr(args, "gff", None)
    scan = validate_gff(Path(gff), drops=drops) if gff else None

    intersects = getattr(args, "intersects", None)
    if intersects:
        paths = [Path(x) for x in intersects]
        ids = _intersect_ids(paths)
        if scan is None:
            print(
                f"[contract] {len(ids)} expression IDs loaded; no --gff given, "
                "so ID recovery could not be checked.",
                file=sys.stderr,
            )
        else:
            check_expression_ids_in_gff(
                ids,
                scan,
                label=f"{len(paths)} intersect file(s)",
                min_overlap=getattr(args, "min_id_overlap", DEFAULT_MIN_ID_OVERLAP),
            )


def _provenance(cmd: str, args: argparse.Namespace, inputs: list[Path]) -> Provenance:
    skip = {"cmd", "no_validate"}
    params = {k: v for k, v in sorted(vars(args).items()) if k not in skip and v is not None}
    return Provenance(command=cmd, inputs=inputs, params=params)


def _dispatch(args: argparse.Namespace, drops: DropLog) -> None:
    if args.cmd == "annotate":
        cis_dir = Path(args.cis_dir)
        prov = _provenance("annotate", args, sorted(cis_dir.glob("*_cis_regulation_module_output.txt")))
        run_annotate(
            cis_dir=cis_dir,
            outdir=Path(args.outdir),
            ncbi_api_key=args.ncbi_api_key,
            delay=args.delay,
            positions=set(args.positions.split(",")),
            use_emapper=args.emapper,
            gff=Path(args.gff) if args.gff else None,
            emapper_db=Path(args.emapper_db) if args.emapper_db else None,
            tax_scope=args.tax_scope,
            emapper_cpu=args.emapper_cpu,
            prov=prov,
        )
        return

    if args.cmd == "exclusive":
        gff_path = Path(args.gff) if args.gff else None
        intersects = [Path(x) for x in args.intersects]
        prov = _provenance("exclusive", args, [Path(args.lnc_ids), *intersects])
        run_exclusive(
            intersect_paths=intersects,
            lncrna_ids_path=Path(args.lnc_ids),
            outdir=Path(args.outdir),
            mode=args.mode,
            min_required=args.min_required,
            level=args.level,
            gff=gff_path,
            prov=prov,
        )
        return

    if args.cmd == "neighbors":
        exclusive_files = [Path(x) for x in args.exclusive]
        intersects = [Path(x) for x in args.intersects]
        prov = _provenance("neighbors", args, [Path(args.gff), *exclusive_files, *intersects])
        run_neighbors(
            gff=Path(args.gff),
            exclusive_files=exclusive_files,
            intersect_paths=intersects,
            out_tsv=Path(args.out),
            k=args.k,
            exclusive_level=args.level,
            prov=prov,
        )
        return

    if args.cmd == "extract-ids":
        out = Path(args.out)
        out.parent.mkdir(parents=True, exist_ok=True)
        prov = _provenance("extract-ids", args, [Path(args.gff)])
        ids, stats = extract_lncrna_ids_from_gff(Path(args.gff), level=args.level)
        with out.open("w", encoding="utf-8") as f:
            prov.write_header(f)
            for x in sorted(ids):
                f.write(x + "\n")
        if args.summary:
            write_extract_ids_summary(Path(args.summary), Path(args.gff), args.level, stats, prov=prov)
        return

    if args.cmd == "intergenic-stats":
        prov = _provenance("intergenic-stats", args, [Path(args.gff)])
        run_intergenic_stats(
            gff=Path(args.gff),
            out_summary=Path(args.summary) if args.summary else None,
            out_distances=Path(args.distances_out) if args.distances_out else None,
            prov=prov,
        )
        return

    if args.cmd == "intron-stats":
        prov = _provenance("intron-stats", args, [Path(args.tsv)])
        run_intron_stats(
            tsv=Path(args.tsv),
            column_index=args.column_index,
            out_summary=Path(args.summary) if args.summary else None,
            prov=prov,
        )
        return

    if args.cmd == "run":
        outdir = Path(args.outdir)
        outdir.mkdir(parents=True, exist_ok=True)
        gff = Path(args.gff)
        intersects = [Path(x) for x in args.intersects]
        prov = _provenance("run", args, [gff, *intersects])

        if args.write_gff_provenance:
            target = write_gff_provenance_file(gff)
            print(f"[provenance] wrote {target}", file=sys.stderr)

        # 1) IDs
        if args.lnc_ids_from_gff:
            suffix = args.level
            lncs_path = outdir / f"lncRNA_IDs.from_gff.{suffix}.txt"
            ids, stats = extract_lncrna_ids_from_gff(gff, level=args.level)
            with lncs_path.open("w", encoding="utf-8") as f:
                prov.write_header(f)
                for x in sorted(ids):
                    f.write(x + "\n")
            write_extract_ids_summary(
                outdir / f"lncRNA_IDs.from_gff.{suffix}.summary.tsv", gff, args.level, stats, prov=prov
            )
        else:
            lncs_path = Path(args.lnc_ids)

        # 2) exclusive
        run_exclusive(
            intersect_paths=intersects,
            lncrna_ids_path=lncs_path,
            outdir=outdir,
            mode=args.mode,
            min_required=args.min_required,
            level=args.level,
            gff=gff if args.level == "gene" else None,
            prov=prov,
        )

        # 3) cis module (neighbors)
        exclusive_files = sorted(outdir.glob("*_exclusive_lncRNAs_*.txt"))
        run_neighbors(
            gff=gff,
            exclusive_files=exclusive_files,
            intersect_paths=intersects,
            out_tsv=outdir / "neighbors.tsv",
            k=args.k,
            exclusive_level=args.level,
            prov=prov,
        )

        run_intergenic_stats(
            gff=gff,
            out_summary=outdir / "intergenic.summary.tsv",
            out_distances=outdir / "intergenic.distances.tsv",
            prov=prov,
        )

        if args.intron_tsv:
            run_intron_stats(
                tsv=Path(args.intron_tsv),
                column_index=args.intron_column_index,
                out_summary=outdir / "intron.summary.tsv",
                prov=prov,
            )
        return

    raise SystemExit(f"unknown command: {args.cmd}")


def main() -> None:
    parser = build_parser()
    args = parser.parse_args()
    drops = DropLog(source=args.cmd)
    set_default_drop_log(drops)

    try:
        _enforce_contract(args, drops)
        _dispatch(args, drops)
    except InputContractError as exc:
        drops.report()
        print(f"ERROR: {exc}", file=sys.stderr)
        raise SystemExit(2)
    except (ValueError, OSError) as exc:
        # The pipeline raises ValueError with a diagnostic message when the
        # inputs are self-consistent but cannot produce a result. That is a
        # user-facing condition, not a crash -- report it without a traceback.
        drops.report()
        print(f"ERROR: {exc}", file=sys.stderr)
        raise SystemExit(1)

    drops.report()


if __name__ == "__main__":
    main()
