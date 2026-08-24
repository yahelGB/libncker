"""The input contract must fail fast, and name the check that failed."""

from __future__ import annotations

from pathlib import Path

import pytest

from lncker.validate import (
    DropLog,
    InputContractError,
    check_expression_ids_in_gff,
    check_gff_version_header,
    read_gff_provenance_pragmas,
    scan_gff,
    validate_gff,
)


def test_valid_fixture_passes_every_check(mini_gff: Path) -> None:
    scan = validate_gff(mini_gff)
    assert scan.drops.total == 0
    assert scan.feature_type_counts["lnc_RNA"] == 3


def test_missing_gff_version_header_is_rejected(fixtures: Path) -> None:
    with pytest.raises(InputContractError) as exc:
        validate_gff(fixtures / "no_version.gff")
    assert exc.value.check == "gff_version_header"


def test_gff_version_2_is_rejected(tmp_path: Path) -> None:
    p = tmp_path / "v2.gff"
    p.write_text("##gff-version 2\nX\ts\tgene\t1\t9\t.\t+\t.\tID=g\n")
    with pytest.raises(InputContractError) as exc:
        check_gff_version_header(p)
    assert exc.value.check == "gff_version_header"


def test_ensembl_style_gff_is_rejected_loudly(fixtures: Path) -> None:
    """The historical silent failure: Ensembl type names yielded empty output."""
    with pytest.raises(InputContractError) as exc:
        validate_gff(fixtures / "ensembl_style.gff")
    assert exc.value.check == "feature_types_present"
    assert "lnc_RNA" in str(exc.value)


def test_dangling_parent_is_rejected(fixtures: Path) -> None:
    with pytest.raises(InputContractError) as exc:
        validate_gff(fixtures / "dangling_parent.gff")
    assert exc.value.check == "parent_resolvable"
    assert "gene-MISSING" in str(exc.value)


def test_missing_file_is_rejected(tmp_path: Path) -> None:
    with pytest.raises(InputContractError) as exc:
        validate_gff(tmp_path / "nope.gff")
    assert exc.value.check == "gff_exists"


def test_expression_id_overlap_is_reported_and_enforced(mini_gff: Path) -> None:
    scan = scan_gff(mini_gff)
    found, total, frac = check_expression_ids_in_gff(
        ["XR_001.1", "XM_003.1"], scan, label="ok", min_overlap=0.95
    )
    assert (found, total, frac) == (2, 2, 1.0)

    with pytest.raises(InputContractError) as exc:
        check_expression_ids_in_gff(
            ["XR_001.1", "NOT_IN_GFF_1", "NOT_IN_GFF_2"], scan, label="bad", min_overlap=0.95
        )
    assert exc.value.check == "expression_ids_in_gff"
    assert "33.33%" in str(exc.value)


def test_overlap_exactly_at_the_floor_passes(mini_gff: Path) -> None:
    scan = scan_gff(mini_gff)
    ids = ["XR_001.1", "XR_002.1", "XR_003.1", "XM_001.1", "NOPE"]
    found, total, frac = check_expression_ids_in_gff(ids, scan, label="edge", min_overlap=0.8)
    assert frac == pytest.approx(0.8)


def test_drop_log_counts_are_exact_beyond_the_warning_cap(capsys) -> None:
    log = DropLog(source="t")
    for i in range(50):
        log.drop("some_reason", f"line {i}")
    assert log.counts["some_reason"] == 50
    assert log.total == 50
    printed = capsys.readouterr().err
    # warnings are capped, the count is not
    assert printed.count("dropped record (some_reason)") == 10
    assert "further 'some_reason' warnings suppressed" in printed


def test_drop_log_summary_rows() -> None:
    log = DropLog(source="t")
    log.drop("a", "x")
    log.drop("b", "y")
    log.drop("a", "z")
    assert log.summary_rows() == [("dropped_total", 3), ("dropped_a", 2), ("dropped_b", 1)]


def test_ncbi_pragmas_are_read(mini_gff: Path) -> None:
    pragmas = read_gff_provenance_pragmas(mini_gff)
    assert pragmas["genome-build"] == "TEST1"
    assert pragmas["genome-build-accession"] == "NCBI_Assembly:GCF_000000000.0"
    assert pragmas["gff-version"] == "3"
