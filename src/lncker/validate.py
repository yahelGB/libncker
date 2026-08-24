"""
Input contract enforcement for lncker.

lncker declares a single supported input format: **NCBI RefSeq GFF3**.
This module validates that contract up front and fails fast with a specific,
named error, instead of letting a foreign annotation silently produce empty
outputs (the historical failure mode -- see CLAUDE.md).

Checks, in the order they run:

    gff_version_header    ``##gff-version 3`` must be the first pragma.
    feature_types_present ``gene``/``mRNA``/``lnc_RNA``/``exon``/``CDS`` all present.
    id_attribute_present  every feature carries ``ID=``.
    parent_resolvable     every ``Parent=`` names a feature that exists.
    expression_ids_in_gff >= ``min_overlap`` of expression-matrix IDs are in the GFF.

Every check that fails raises :class:`InputContractError` naming the check.
Malformed records are never skipped silently: they go through :class:`DropLog`,
which WARNs and counts, and whose total is printed at the end of the run.
"""

from __future__ import annotations

import sys
from collections import Counter
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, Iterable, Iterator, List, Optional, Set, Tuple

from lncker.io_utils import open_text

# Feature types an NCBI RefSeq GFF3 of an annotated genome must contain.
REQUIRED_FEATURE_TYPES = ("gene", "mRNA", "lnc_RNA", "exon", "CDS")

# Transcript-level feature types whose IDs may appear in an expression matrix.
TRANSCRIPT_FEATURE_TYPES = (
    "mRNA",
    "lnc_RNA",
    "transcript",
    "tRNA",
    "rRNA",
    "snRNA",
    "snoRNA",
    "guide_RNA",
    "primary_transcript",
    "antisense_RNA",
    "V_gene_segment",
    "C_gene_segment",
)

# Per-reason cap on individual WARN lines. Counts stay exact regardless; only
# the printing is capped, so a systematically broken file cannot flood stderr.
MAX_WARNINGS_PER_REASON = 10

DEFAULT_MIN_ID_OVERLAP = 0.95


class InputContractError(ValueError):
    """Raised when an input violates lncker's declared input contract."""

    def __init__(self, check: str, detail: str) -> None:
        self.check = check
        self.detail = detail
        super().__init__(f"input contract check '{check}' FAILED: {detail}")


@dataclass
class DropLog:
    """Counts and reports every record lncker refuses to use.

    ``lncker`` must never drop a record without saying so. Callers report a
    drop with :meth:`drop`; :meth:`report` prints the per-reason table and the
    grand total, and :meth:`summary_rows` feeds the same numbers into the
    command's ``summary.tsv``.
    """

    source: str
    counts: Counter = field(default_factory=Counter)
    _warned: Counter = field(default_factory=Counter)
    _seen: Set[Tuple[str, str]] = field(default_factory=set)

    def drop(self, reason: str, detail: str) -> None:
        # lncker reads a GFF in several passes; the same bad record must be
        # counted once, not once per pass.
        key = (reason, detail)
        if key in self._seen:
            return
        self._seen.add(key)
        self.counts[reason] += 1
        if self._warned[reason] < MAX_WARNINGS_PER_REASON:
            self._warned[reason] += 1
            print(f"WARN [{self.source}] dropped record ({reason}): {detail}", file=sys.stderr)
            if self._warned[reason] == MAX_WARNINGS_PER_REASON:
                print(
                    f"WARN [{self.source}] further '{reason}' warnings suppressed; "
                    f"the final count is exact.",
                    file=sys.stderr,
                )

    @property
    def total(self) -> int:
        return sum(self.counts.values())

    def report(self) -> None:
        if not self.total:
            print(f"[{self.source}] records dropped: 0", file=sys.stderr)
            return
        print(f"[{self.source}] records dropped: {self.total}", file=sys.stderr)
        for reason, n in sorted(self.counts.items(), key=lambda kv: -kv[1]):
            print(f"[{self.source}]   {reason}\t{n}", file=sys.stderr)

    def summary_rows(self, prefix: str = "dropped") -> List[Tuple[str, int]]:
        rows = [(f"{prefix}_total", self.total)]
        rows.extend((f"{prefix}_{reason}", n) for reason, n in sorted(self.counts.items()))
        return rows


_DEFAULT_DROPS: Optional[DropLog] = None


def set_default_drop_log(log: DropLog) -> None:
    """Route drops from every parser in this process into one log."""
    global _DEFAULT_DROPS
    _DEFAULT_DROPS = log


def default_drop_log() -> DropLog:
    global _DEFAULT_DROPS
    if _DEFAULT_DROPS is None:
        _DEFAULT_DROPS = DropLog(source="lncker")
    return _DEFAULT_DROPS


def parse_gff_record(
    line: str, lineno: int, drops: DropLog
) -> Optional[Tuple[str, str, str, int, int, str, str, str]]:
    """Parse one GFF3 data line, or report it to ``drops`` and return None.

    This is the only place lncker turns a GFF line into fields, so a malformed
    record is dropped identically -- and counted -- no matter which command is
    reading the file.
    """
    parts = line.rstrip("\n").split("\t")
    if len(parts) < 9:
        drops.drop("malformed_column_count", f"line {lineno}: {len(parts)} columns, expected 9")
        return None
    contig, src, ftype, start_s, end_s, score, strand, phase, attrs_s = parts[:9]
    try:
        start = int(start_s)
        end = int(end_s)
    except ValueError:
        drops.drop("non_integer_coordinates", f"line {lineno}: start={start_s!r} end={end_s!r}")
        return None
    if start < 1 or end < 1:
        drops.drop("non_positive_coordinates", f"line {lineno}: start={start} end={end}")
        return None
    if start > end:
        drops.drop("start_after_end", f"line {lineno}: start={start} > end={end}")
        return None
    return contig, src, ftype, start, end, strand, phase, attrs_s


def _split_parents(raw: str) -> List[str]:
    """GFF3 ``Parent`` is a comma-separated list; NCBI uses a single value."""
    return [p for p in (x.strip() for x in raw.split(",")) if p]


def check_gff_version_header(gff: Path) -> str:
    """``##gff-version 3`` must be the first non-blank line."""
    with open_text(gff) as f:
        for line in f:
            if not line.strip():
                continue
            first = line.strip()
            break
        else:
            raise InputContractError("gff_version_header", f"{gff} is empty")

    if not first.startswith("##gff-version"):
        raise InputContractError(
            "gff_version_header",
            f"{gff} does not start with '##gff-version'; first line was {first!r}. "
            "lncker requires an NCBI RefSeq GFF3.",
        )
    version = first.split(None, 1)[1].strip() if len(first.split(None, 1)) > 1 else ""
    if not version.startswith("3"):
        raise InputContractError(
            "gff_version_header",
            f"{gff} declares GFF version {version!r}; lncker requires GFF3.",
        )
    return first


def read_gff_provenance_pragmas(gff: Path) -> Dict[str, str]:
    """Collect the ``#!key value`` pragmas NCBI writes into its GFF headers."""
    out: Dict[str, str] = {}
    with open_text(gff) as f:
        for line in f:
            if not line.startswith("#"):
                break
            s = line.rstrip("\n")
            if s.startswith("#!"):
                parts = s[2:].split(None, 1)
                if len(parts) == 2:
                    out[parts[0]] = parts[1].strip()
            elif s.startswith("##gff-version"):
                out["gff-version"] = s.split(None, 1)[1].strip()
    return out


@dataclass
class GffScan:
    """Everything a single pass over the GFF can tell us."""

    feature_type_counts: Counter
    total_records: int
    ids: Set[str]
    transcript_ids: Set[str]
    unresolved_parents: Counter
    mrna_per_gene: Counter
    contigs: Set[str]
    genes_per_contig: Counter
    features_per_contig: Counter
    drops: DropLog


def scan_gff(gff: Path, drops: Optional[DropLog] = None) -> GffScan:
    """One strict pass over the GFF, collecting counts and integrity data.

    Unlike :func:`lncker.gff_utils.iter_gff`, a record that cannot be parsed is
    reported through the drop log rather than skipped in silence.
    """
    drops = drops if drops is not None else default_drop_log()

    ftypes: Counter = Counter()
    ids: Set[str] = set()
    transcript_ids: Set[str] = set()
    parents: Counter = Counter()
    mrna_per_gene: Counter = Counter()
    contigs: Set[str] = set()
    genes_per_contig: Counter = Counter()
    features_per_contig: Counter = Counter()
    total = 0

    with open_text(gff) as f:
        for lineno, line in enumerate(f, start=1):
            if line.startswith("#") or not line.strip():
                continue
            rec = parse_gff_record(line, lineno, drops)
            if rec is None:
                continue
            contig, _src, ftype, start, end, _strand, _phase, attrs_s = rec

            attrs: Dict[str, str] = {}
            for item in attrs_s.split(";"):
                if "=" in item:
                    k, v = item.split("=", 1)
                    attrs[k] = v

            total += 1
            ftypes[ftype] += 1
            contigs.add(contig)
            features_per_contig[contig] += 1

            fid = attrs.get("ID", "")
            if fid:
                ids.add(fid)
            if ftype in TRANSCRIPT_FEATURE_TYPES and fid:
                transcript_ids.add(fid[4:] if fid.startswith("rna-") else fid)
            if ftype == "gene":
                genes_per_contig[contig] += 1
            for parent in _split_parents(attrs.get("Parent", "")):
                parents[parent] += 1
                if ftype == "mRNA":
                    mrna_per_gene[parent] += 1

    unresolved = Counter({p: n for p, n in parents.items() if p not in ids})

    return GffScan(
        feature_type_counts=ftypes,
        total_records=total,
        ids=ids,
        transcript_ids=transcript_ids,
        unresolved_parents=unresolved,
        mrna_per_gene=mrna_per_gene,
        contigs=contigs,
        genes_per_contig=genes_per_contig,
        features_per_contig=features_per_contig,
        drops=drops,
    )


def check_feature_types_present(scan: GffScan, gff: Path) -> None:
    missing = [t for t in REQUIRED_FEATURE_TYPES if not scan.feature_type_counts.get(t)]
    if missing:
        seen = ", ".join(f"{t}={n}" for t, n in scan.feature_type_counts.most_common(10))
        raise InputContractError(
            "feature_types_present",
            f"{gff} has no {missing} feature(s). lncker requires NCBI RefSeq type names; "
            f"Ensembl-style GFFs use 'transcript'/'ncRNA' and are not supported. "
            f"Types found: {seen or '(none)'}",
        )


def check_id_attribute_present(scan: GffScan, gff: Path) -> None:
    if not scan.ids:
        raise InputContractError(
            "id_attribute_present",
            f"{gff} has no feature carrying an 'ID=' attribute.",
        )


def check_parent_resolvable(scan: GffScan, gff: Path) -> None:
    if scan.unresolved_parents:
        worst = ", ".join(f"{p} (x{n})" for p, n in scan.unresolved_parents.most_common(5))
        raise InputContractError(
            "parent_resolvable",
            f"{gff} has {len(scan.unresolved_parents)} Parent value(s) that name no feature ID. "
            f"Examples: {worst}",
        )


def check_expression_ids_in_gff(
    expression_ids: Iterable[str],
    scan: GffScan,
    label: str,
    min_overlap: float = DEFAULT_MIN_ID_OVERLAP,
) -> Tuple[int, int, float]:
    """Report what fraction of expression-matrix IDs exist in the GFF; abort below the floor."""
    ids = {x for x in (s.strip() for s in expression_ids) if x}
    if not ids:
        raise InputContractError("expression_ids_in_gff", f"{label} contains no usable IDs.")

    found = len(ids & scan.transcript_ids)
    frac = found / len(ids)
    print(
        f"[contract] {label}: {found}/{len(ids)} IDs found in GFF ({100 * frac:.2f}%)",
        file=sys.stderr,
    )
    if frac < min_overlap:
        sample = sorted(ids - scan.transcript_ids)[:5]
        raise InputContractError(
            "expression_ids_in_gff",
            f"{label}: only {100 * frac:.2f}% of IDs are present in the GFF "
            f"(required: {100 * min_overlap:.2f}%). This normally means the expression matrix "
            f"and the GFF come from different annotation releases, or the IDs are at the wrong "
            f"level (gene vs transcript). Missing examples: {sample}",
        )
    return found, len(ids), frac


def validate_gff(gff: Path, drops: Optional[DropLog] = None) -> GffScan:
    """Run every GFF-only contract check and return the scan for reuse."""
    if not gff.exists():
        raise InputContractError("gff_exists", f"{gff} does not exist")
    check_gff_version_header(gff)
    scan = scan_gff(gff, drops=drops)
    check_feature_types_present(scan, gff)
    check_id_attribute_present(scan, gff)
    check_parent_resolvable(scan, gff)
    print(
        f"[contract] {gff.name}: OK "
        f"({scan.total_records} records, {len(scan.contigs)} sequences, "
        f"{scan.feature_type_counts.get('gene', 0)} genes)",
        file=sys.stderr,
    )
    return scan
