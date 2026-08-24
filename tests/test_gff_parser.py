"""Parser conformance: what lncker parses must equal what is in the file."""

from __future__ import annotations

import collections
import subprocess
from pathlib import Path

import pytest

from lncker.gff_utils import (
    build_lnc_gene_to_transcripts,
    build_transcript_to_gene_map,
    extract_lncrna_ids_from_gff,
    iter_gff,
    norm_rna_id,
    parse_attrs,
)
from lncker.validate import DropLog, scan_gff


def _grep_counts(gff: Path) -> dict[str, int]:
    """Ground truth straight from the file, computed without lncker."""
    counts: collections.Counter = collections.Counter()
    for line in gff.read_text().splitlines():
        if line.startswith("#") or not line.strip():
            continue
        parts = line.split("\t")
        if len(parts) < 9:
            continue
        try:
            start, end = int(parts[3]), int(parts[4])
        except ValueError:
            continue
        if start < 1 or end < 1 or start > end:
            continue
        counts[parts[2]] += 1
    return dict(counts)


def test_parsed_feature_counts_match_the_raw_file(mini_gff: Path) -> None:
    parsed = collections.Counter(f.ftype for f in iter_gff(mini_gff))
    assert dict(parsed) == _grep_counts(mini_gff)


def test_scan_counts_match_the_raw_file(mini_gff: Path) -> None:
    scan = scan_gff(mini_gff)
    assert dict(scan.feature_type_counts) == _grep_counts(mini_gff)
    assert scan.total_records == sum(_grep_counts(mini_gff).values())
    assert scan.drops.total == 0


def test_scaffold_level_counts(mini_gff: Path) -> None:
    scan = scan_gff(mini_gff)
    assert len(scan.contigs) == 3
    assert scan.genes_per_contig["NW_TEST01.1"] == 4
    assert scan.genes_per_contig["NW_TEST02.1"] == 1
    assert scan.genes_per_contig["NW_TEST03.1"] == 2
    # NW_TEST02.1 is the single-gene scaffold the fixture exists to cover.
    single = [c for c, n in scan.genes_per_contig.items() if n < 2]
    assert single == ["NW_TEST02.1"]


def test_multi_isoform_gene_is_kept_as_two_transcripts(mini_gff: Path) -> None:
    scan = scan_gff(mini_gff)
    assert scan.mrna_per_gene["gene-LOC001"] == 2
    tx2gene = build_transcript_to_gene_map(mini_gff)
    assert tx2gene["XM_001.1"] == "gene-LOC001"
    assert tx2gene["XM_002.1"] == "gene-LOC001"


def test_features_at_scaffold_edges_are_parsed(mini_gff: Path) -> None:
    feats = [f for f in iter_gff(mini_gff) if f.contig == "NW_TEST03.1"]
    starts = {f.start for f in feats}
    ends = {f.end for f in feats}
    assert 1 in starts, "feature starting at base 1 was dropped"
    assert 1000 in ends, "feature ending at the last base was dropped"


def test_overlapping_antisense_features_are_both_kept(mini_gff: Path) -> None:
    by_id = {f.attrs.get("ID"): f for f in iter_gff(mini_gff)}
    mrna = by_id["rna-XM_003.1"]
    lnc = by_id["rna-XR_002.1"]
    assert mrna.strand == "+" and lnc.strand == "-"
    assert lnc.start <= mrna.end and mrna.start <= lnc.end, "fixture must overlap"


def test_lncrna_discovery_at_both_levels(mini_gff: Path) -> None:
    tx, stats = extract_lncrna_ids_from_gff(mini_gff, level="transcript")
    assert tx == {"XR_001.1", "XR_002.1", "XR_003.1"}
    genes, _ = extract_lncrna_ids_from_gff(mini_gff, level="gene")
    assert genes == {"gene-LOC002", "gene-LOC004", "gene-LOC006"}
    assert stats["transcripts_lncRNA"] == 3


def test_lncrna_discovery_works_without_gene_biotype(tmp_path: Path, mini_gff: Path) -> None:
    """Discovery unions gene_biotype with the Parent of every lnc_RNA."""
    stripped = tmp_path / "no_biotype.gff"
    stripped.write_text(mini_gff.read_text().replace(";gene_biotype=lncRNA", ""))
    genes, _ = extract_lncrna_ids_from_gff(stripped, level="gene")
    assert genes == {"gene-LOC002", "gene-LOC004", "gene-LOC006"}


def test_gene_to_transcript_expansion(mini_gff: Path) -> None:
    mapping = build_lnc_gene_to_transcripts(mini_gff)
    assert mapping["gene-LOC002"] == {"XR_001.1"}
    assert mapping["gene-LOC006"] == {"XR_003.1"}


def test_norm_rna_id_strips_only_the_rna_prefix() -> None:
    assert norm_rna_id("rna-XR_001.1") == "XR_001.1"
    assert norm_rna_id("gene-LOC001") == "gene-LOC001"
    assert norm_rna_id("XR_001.1") == "XR_001.1"


def test_parse_attrs_keeps_values_containing_equals() -> None:
    got = parse_attrs("ID=a;note=x=y;gbkey=Gene")
    assert got == {"ID": "a", "note": "x=y", "gbkey": "Gene"}


def test_malformed_records_are_counted_not_silently_skipped(fixtures: Path) -> None:
    drops = DropLog(source="test")
    scan = scan_gff(fixtures / "malformed.gff", drops=drops)
    assert scan.drops.total == 4
    assert scan.drops.counts["non_integer_coordinates"] == 1
    assert scan.drops.counts["start_after_end"] == 1
    assert scan.drops.counts["non_positive_coordinates"] == 1
    assert scan.drops.counts["malformed_column_count"] == 1
    # the six well-formed records survive: region, gene, mRNA, exon, CDS, lnc_RNA
    assert scan.total_records == 6
    assert dict(scan.feature_type_counts) == {
        "region": 1, "gene": 1, "mRNA": 1, "exon": 1, "CDS": 1, "lnc_RNA": 1,
    }
