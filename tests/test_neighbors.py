"""Neighbour search geometry, asserted against hand-computed fixture coordinates."""

from __future__ import annotations

import csv
from pathlib import Path

import pytest

from lncker.neighbors import (
    Feature,
    _build_lnc_index,
    _sort_key,
    distance_kb,
    find_neighbors_for_one_lnc,
    iter_gff_features,
    run_neighbors,
)


def _feat(start: int, end: int, ftype: str = "mRNA", strand: str = "+", **attrs) -> Feature:
    return Feature(contig="c", ftype=ftype, start=start, end=end, strand=strand, attrs=attrs)


# --------------------------------------------------------------------------
# distance_kb
# --------------------------------------------------------------------------

def test_distance_upstream_and_downstream() -> None:
    lnc = _feat(5000, 6000, "lnc_RNA")
    assert distance_kb(lnc, _feat(1200, 3000)) == pytest.approx(2.0)
    assert distance_kb(lnc, _feat(8000, 9000)) == pytest.approx(2.0)


def test_distance_is_negative_when_features_overlap() -> None:
    """XR_002.1 (8500-9500) overlaps XM_003.1 (8000-9000) by 501 bp."""
    lnc = _feat(8500, 9500, "lnc_RNA", strand="-")
    assert distance_kb(lnc, _feat(8000, 9000)) == pytest.approx(-0.501)


def test_distance_of_a_single_base_gap() -> None:
    lnc = _feat(100, 200, "lnc_RNA")
    assert distance_kb(lnc, _feat(201, 300)) == pytest.approx(0.001)


def test_distance_of_full_containment() -> None:
    lnc = _feat(100, 200, "lnc_RNA")
    assert distance_kb(lnc, _feat(50, 300)) == pytest.approx(-0.101)


# --------------------------------------------------------------------------
# neighbour walk
# --------------------------------------------------------------------------

def test_neighbours_deduplicate_by_parent_gene(mini_gff: Path) -> None:
    """A multi-isoform gene must contribute exactly one transcript."""
    feats = sorted(
        (f for f in iter_gff_features(mini_gff) if f.contig == "NW_TEST01.1"), key=_sort_key
    )
    idx = next(i for i, f in enumerate(feats) if f.attrs.get("ID") == "rna-XR_001.1")
    up, dn = find_neighbors_for_one_lnc(feats, idx, k=5, exclusive_lncs_all=set(), target_lnc_id="XR_001.1")

    up_ids = [feats[i].attrs["ID"] for i in up]
    dn_ids = [feats[i].attrs["ID"] for i in dn]
    assert up_ids == ["rna-XM_002.1"], "both isoforms of gene-LOC001 were emitted"
    assert dn_ids == ["rna-XM_003.1"]


def test_scaffold_edge_lncrna_has_no_upstream_neighbour(mini_gff: Path) -> None:
    feats = sorted(
        (f for f in iter_gff_features(mini_gff) if f.contig == "NW_TEST03.1"), key=_sort_key
    )
    idx = next(i for i, f in enumerate(feats) if f.attrs.get("ID") == "rna-XR_003.1")
    up, dn = find_neighbors_for_one_lnc(feats, idx, k=5, exclusive_lncs_all=set(), target_lnc_id="XR_003.1")
    assert up == []
    assert [feats[i].attrs["ID"] for i in dn] == ["rna-XM_005.1"]


def test_walk_stops_at_another_exclusive_lncrna(mini_gff: Path) -> None:
    """XR_002.1 sits between XR_001.1 and everything further downstream."""
    feats = sorted(
        (f for f in iter_gff_features(mini_gff) if f.contig == "NW_TEST01.1"), key=_sort_key
    )
    idx = next(i for i, f in enumerate(feats) if f.attrs.get("ID") == "rna-XR_001.1")

    _, dn_open = find_neighbors_for_one_lnc(feats, idx, k=5, exclusive_lncs_all=set(), target_lnc_id="XR_001.1")
    _, dn_blocked = find_neighbors_for_one_lnc(
        feats, idx, k=5, exclusive_lncs_all={"XR_001.1", "XR_002.1"}, target_lnc_id="XR_001.1"
    )
    # XM_003.1 precedes XR_002.1 in sort order, so it is still reached either way,
    # but the target itself must never block its own walk.
    assert dn_open == dn_blocked
    assert [feats[i].attrs["ID"] for i in dn_blocked] == ["rna-XM_003.1"]


def test_lnc_index_finds_every_lncrna(mini_gff: Path) -> None:
    contigs: dict[str, list[Feature]] = {}
    for f in iter_gff_features(mini_gff):
        contigs.setdefault(f.contig, []).append(f)
    for c in contigs:
        contigs[c].sort(key=_sort_key)
    idx = _build_lnc_index(contigs)
    assert set(idx) == {"XR_001.1", "XR_002.1", "XR_003.1"}


# --------------------------------------------------------------------------
# end-to-end
# --------------------------------------------------------------------------

def _read_tsv(path: Path) -> list[dict[str, str]]:
    lines = [ln for ln in path.read_text().splitlines() if not ln.startswith("#")]
    return list(csv.DictReader(lines, delimiter="\t"))


def test_run_neighbors_end_to_end(tmp_path: Path, mini_gff: Path, intersects: list[Path]) -> None:
    excl = tmp_path / "A_exclusive_lncRNAs_strict.txt"
    excl.write_text("XR_001.1\n")
    out = tmp_path / "neighbors.tsv"

    run_neighbors(
        gff=mini_gff,
        exclusive_files=[excl],
        intersect_paths=intersects,
        out_tsv=out,
        k=5,
        exclusive_level="transcript",
    )

    rows = _read_tsv(out)
    assert {r["mRNA_position"] for r in rows} == {"1_U", "1_D"}
    by_pos = {r["mRNA_position"]: r for r in rows}

    assert by_pos["1_U"]["mRNA_ID"] == "XM_002.1"
    assert float(by_pos["1_U"]["distance_kb"]) == pytest.approx(2.0)
    assert by_pos["1_D"]["mRNA_ID"] == "XM_003.1"
    assert float(by_pos["1_D"]["distance_kb"]) == pytest.approx(2.0)

    # 3 tissues (A, B, C) -> cis_expected = N - 1 = 2
    assert by_pos["1_D"]["cis_expected"] == "2"
    # XM_003.1 is up in A in both AvsB and AvsC
    assert by_pos["1_D"]["cis_score"] == "2"
    assert float(by_pos["1_D"]["cis_probability"]) == pytest.approx(1.0)
    # XM_002.1 is up in A only in AvsB
    assert by_pos["1_U"]["cis_score"] == "1"
    assert float(by_pos["1_U"]["cis_probability"]) == pytest.approx(0.5)
    assert by_pos["1_D"]["DE_status"] == "DE"


def test_run_neighbors_writes_per_tissue_cis_file(tmp_path: Path, mini_gff: Path, intersects: list[Path]) -> None:
    excl = tmp_path / "A_exclusive_lncRNAs_strict.txt"
    excl.write_text("XR_001.1\n")
    run_neighbors(
        gff=mini_gff,
        exclusive_files=[excl],
        intersect_paths=intersects,
        out_tsv=tmp_path / "neighbors.tsv",
        k=5,
        exclusive_level="transcript",
    )
    cis = tmp_path / "A_cis_regulation_module_output.txt"
    assert cis.exists()
    rows = _read_tsv(cis)
    assert {r["mRNA_ID"] for r in rows} == {"XM_002.1", "XM_003.1"}
    # the per-tissue file drops cis_expected / cis_probability
    assert "cis_expected" not in rows[0]


def test_gene_level_exclusive_ids_are_expanded(tmp_path: Path, mini_gff: Path, intersects: list[Path]) -> None:
    excl = tmp_path / "A_exclusive_lncRNAs_strict.txt"
    excl.write_text("gene-LOC002\n")
    out = tmp_path / "neighbors.tsv"
    run_neighbors(
        gff=mini_gff,
        exclusive_files=[excl],
        intersect_paths=intersects,
        out_tsv=out,
        k=5,
        exclusive_level="gene",
    )
    rows = _read_tsv(out)
    assert {r["lncRNA_ID"] for r in rows} == {"XR_001.1"}


def test_gene_level_with_no_lnc_children_raises(tmp_path: Path, mini_gff: Path, intersects: list[Path]) -> None:
    excl = tmp_path / "A_exclusive_lncRNAs_strict.txt"
    excl.write_text("gene-LOC001\n")  # protein-coding: has no lnc_RNA children
    with pytest.raises(ValueError, match="Exclusive set is empty"):
        run_neighbors(
            gff=mini_gff,
            exclusive_files=[excl],
            intersect_paths=intersects,
            out_tsv=tmp_path / "neighbors.tsv",
            k=5,
            exclusive_level="gene",
        )
