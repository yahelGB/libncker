"""Tissue-exclusive selection, strict and lenient."""

from __future__ import annotations

from pathlib import Path

import pytest

from lncker.exclusive import run_exclusive
from lncker.expr import parse_expression


def _ids(path: Path) -> set[str]:
    return {
        ln.strip()
        for ln in path.read_text().splitlines()
        if ln.strip() and not ln.startswith("#")
    }


def _summary(outdir: Path) -> dict[str, str]:
    out = {}
    for ln in (outdir / "summary.tsv").read_text().splitlines():
        if ln.startswith("#") or "\t" not in ln:
            continue
        k, v = ln.split("\t", 1)
        out[k] = v
    return out


def test_parse_expression_both_orientations() -> None:
    assert parse_expression("Up_A_Down_B") == parse_expression("Down_B_Up_A")
    p = parse_expression("Up_A_Down_B")
    assert (p.up, p.down) == ("A", "B")


def test_parse_expression_tolerates_underscores_in_tissue_names() -> None:
    p = parse_expression("Up_hepato_pancreas_Down_gill_tissue")
    assert (p.up, p.down) == ("hepato_pancreas", "gill_tissue")


@pytest.mark.parametrize("bad", ["", "garbage", "Up_A", "Down_B", "Up__Down_B"])
def test_parse_expression_rejects_malformed(bad: str) -> None:
    assert parse_expression(bad) is None


def test_strict_mode_selects_only_fully_consistent_lncrnas(
    tmp_path: Path, fixtures: Path, intersects: list[Path]
) -> None:
    lnc_ids = tmp_path / "lnc.txt"
    lnc_ids.write_text("XR_001.1\nXR_002.1\nXR_003.1\n")

    run_exclusive(
        intersect_paths=intersects,
        lncrna_ids_path=lnc_ids,
        outdir=tmp_path,
        mode="strict",
        level="transcript",
    )

    # XR_001.1 is up in A in both comparisons involving A -> exclusive
    assert _ids(tmp_path / "A_exclusive_lncRNAs_strict.txt") == {"XR_001.1"}
    # XR_002.1 is up in A once only -> consistent but incomplete
    assert _ids(tmp_path / "A_up_consistent_but_incomplete.txt") == {"XR_002.1"}
    # XR_003.1 is up in B and in C -> mixed direction, excluded everywhere
    assert _ids(tmp_path / "B_exclusive_lncRNAs_strict.txt") == set()
    assert _ids(tmp_path / "C_exclusive_lncRNAs_strict.txt") == set()

    s = _summary(tmp_path)
    assert s["tissues_inferred"] == "3"
    assert s["required_patterns_strict"] == "2"
    assert s["ids_mixed_up_direction"] == "1"
    assert s["exclusive_count_A"] == "1"


def test_lenient_mode_admits_the_near_miss(tmp_path: Path, intersects: list[Path]) -> None:
    lnc_ids = tmp_path / "lnc.txt"
    lnc_ids.write_text("XR_001.1\nXR_002.1\nXR_003.1\n")

    run_exclusive(
        intersect_paths=intersects,
        lncrna_ids_path=lnc_ids,
        outdir=tmp_path,
        mode="lenient",
        min_required=1,
        level="transcript",
    )
    assert _ids(tmp_path / "A_exclusive_lncRNAs_lenient.txt") == {"XR_001.1", "XR_002.1"}
    # lenient mode writes no near-miss file
    assert not (tmp_path / "A_up_consistent_but_incomplete.txt").exists()


def test_summary_counts_are_not_thrown_off_by_provenance_headers(
    tmp_path: Path, intersects: list[Path]
) -> None:
    """Regression: the count used to come from wc -l of the output file."""
    from lncker.provenance import Provenance

    lnc_ids = tmp_path / "lnc.txt"
    lnc_ids.write_text("XR_001.1\nXR_002.1\nXR_003.1\n")
    run_exclusive(
        intersect_paths=intersects,
        lncrna_ids_path=lnc_ids,
        outdir=tmp_path,
        mode="strict",
        level="transcript",
        prov=Provenance(command="test", params={"k": 1}),
    )
    out = tmp_path / "A_exclusive_lncRNAs_strict.txt"
    assert out.read_text().startswith("# lncker provenance")
    assert _summary(tmp_path)["exclusive_count_A"] == "1"
    assert _ids(out) == {"XR_001.1"}


def test_mismatched_id_level_raises_a_diagnostic_error(tmp_path: Path, intersects: list[Path]) -> None:
    lnc_ids = tmp_path / "lnc.txt"
    lnc_ids.write_text("gene-LOC002\ngene-LOC004\n")  # gene IDs, transcript-level run
    with pytest.raises(ValueError, match="No rows from intersect files matched"):
        run_exclusive(
            intersect_paths=intersects,
            lncrna_ids_path=lnc_ids,
            outdir=tmp_path,
            mode="strict",
            level="transcript",
        )


def test_gene_level_without_gff_is_rejected(tmp_path: Path, intersects: list[Path]) -> None:
    lnc_ids = tmp_path / "lnc.txt"
    lnc_ids.write_text("gene-LOC002\n")
    with pytest.raises(ValueError, match="you must provide --gff"):
        run_exclusive(
            intersect_paths=intersects,
            lncrna_ids_path=lnc_ids,
            outdir=tmp_path,
            mode="strict",
            level="gene",
        )


def test_missing_intersect_header_fails_loudly(tmp_path: Path) -> None:
    bad = tmp_path / "bad_intersect.txt"
    bad.write_text("name\tvalue\nXR_001.1\tUp_A_Down_B\n")
    lnc_ids = tmp_path / "lnc.txt"
    lnc_ids.write_text("XR_001.1\n")
    with pytest.raises(ValueError, match="missing ID/Expression header"):
        run_exclusive(
            intersect_paths=[bad],
            lncrna_ids_path=lnc_ids,
            outdir=tmp_path,
            mode="strict",
            level="transcript",
        )
