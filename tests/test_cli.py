"""Every subcommand: --help, plus one minimal real invocation."""

from __future__ import annotations

from pathlib import Path

import pytest

SUBCOMMANDS = [
    "exclusive",
    "neighbors",
    "extract-ids",
    "intergenic-stats",
    "intron-stats",
    "annotate",
    "run",
]


@pytest.mark.parametrize("cmd", SUBCOMMANDS)
def test_help_exits_zero(cli, cmd: str) -> None:
    r = cli(cmd, "--help")
    assert r.returncode == 0, r.stderr
    assert cmd.split("-")[0] in r.stdout or "usage:" in r.stdout


def test_bare_invocation_requires_a_subcommand(cli) -> None:
    r = cli()
    assert r.returncode != 0
    assert "arguments are required: cmd" in r.stderr


def test_unknown_subcommand_is_rejected(cli) -> None:
    r = cli("nope")
    assert r.returncode != 0


# --------------------------------------------------------------------------
# minimal real invocations
# --------------------------------------------------------------------------

def test_extract_ids(cli, tmp_path: Path, mini_gff: Path) -> None:
    out = tmp_path / "ids.txt"
    r = cli("extract-ids", "--gff", str(mini_gff), "--level", "transcript",
            "--out", str(out), "--summary", str(tmp_path / "ids.summary.tsv"))
    assert r.returncode == 0, r.stderr
    ids = {ln for ln in out.read_text().splitlines() if not ln.startswith("#")}
    assert ids == {"XR_001.1", "XR_002.1", "XR_003.1"}
    assert out.read_text().startswith("# lncker provenance")


def test_exclusive(cli, tmp_path: Path, fixtures: Path) -> None:
    lnc = tmp_path / "lnc.txt"
    lnc.write_text("XR_001.1\nXR_002.1\nXR_003.1\n")
    r = cli("exclusive", "--lnc-ids", str(lnc), "--outdir", str(tmp_path),
            "--mode", "strict", "--level", "transcript",
            "--intersects", *[str(p) for p in sorted(fixtures.glob("*vs*_intersect.txt"))])
    assert r.returncode == 0, r.stderr
    assert (tmp_path / "A_exclusive_lncRNAs_strict.txt").exists()
    assert (tmp_path / "summary.tsv").exists()


def test_neighbors(cli, tmp_path: Path, mini_gff: Path, fixtures: Path) -> None:
    excl = tmp_path / "A_exclusive_lncRNAs_strict.txt"
    excl.write_text("XR_001.1\n")
    out = tmp_path / "neighbors.tsv"
    r = cli("neighbors", "--gff", str(mini_gff), "--exclusive", str(excl),
            "--out", str(out), "--k", "5", "--level", "transcript",
            "--intersects", *[str(p) for p in sorted(fixtures.glob("*vs*_intersect.txt"))])
    assert r.returncode == 0, r.stderr
    assert out.exists()
    assert (tmp_path / "neighbors.summary.tsv").exists()


def test_intergenic_stats(cli, tmp_path: Path, mini_gff: Path) -> None:
    r = cli("intergenic-stats", "--gff", str(mini_gff),
            "--summary", str(tmp_path / "ig.tsv"),
            "--distances-out", str(tmp_path / "d.tsv"))
    assert r.returncode == 0, r.stderr
    # NW_TEST01.1 contributes 2 gaps: LOC003/LOC004 overlap (antisense) and are
    # excluded. NW_TEST02.1 is a single-gene scaffold and contributes 0.
    # NW_TEST03.1 contributes 1.
    body = [ln for ln in (tmp_path / "d.tsv").read_text().splitlines() if not ln.startswith("#")]
    assert len(body) - 1 == 3
    assert "gene-LOC003\tgene-LOC004" not in "\n".join(body), "overlapping antisense pair leaked in"


def test_intron_stats(cli, tmp_path: Path) -> None:
    tsv = tmp_path / "introns.tsv"
    tsv.write_text("a\tb\n100\t200\n300\t400\n")
    r = cli("intron-stats", "--tsv", str(tsv), "--column-index", "1",
            "--summary", str(tmp_path / "intron.tsv"))
    assert r.returncode == 0, r.stderr
    assert "Geometric mean" in r.stdout


def test_run_full_pipeline(cli, tmp_path: Path, mini_gff: Path, fixtures: Path) -> None:
    outdir = tmp_path / "out"
    r = cli("run", "--gff", str(mini_gff), "--lnc-ids-from-gff", "--level", "transcript",
            "--outdir", str(outdir), "--mode", "strict", "--k", "5",
            "--intersects", *[str(p) for p in sorted(fixtures.glob("*vs*_intersect.txt"))])
    assert r.returncode == 0, r.stderr
    for name in ("neighbors.tsv", "neighbors.summary.tsv", "summary.tsv",
                 "intergenic.summary.tsv", "lncRNA_IDs.from_gff.transcript.txt"):
        assert (outdir / name).exists(), f"{name} missing"


def test_annotate_reports_an_empty_input_directory(cli, tmp_path: Path) -> None:
    """A real invocation that needs no network: annotate with nothing to do."""
    r = cli("annotate", "--cis-dir", str(tmp_path), "--outdir", str(tmp_path / "ann"))
    assert r.returncode != 0
    assert "No *_cis_regulation_module_output.txt files found" in (r.stderr + r.stdout)


# --------------------------------------------------------------------------
# input contract at the CLI boundary
# --------------------------------------------------------------------------

def test_ensembl_gff_aborts_with_exit_code_2(cli, tmp_path: Path, fixtures: Path) -> None:
    r = cli("extract-ids", "--gff", str(fixtures / "ensembl_style.gff"), "--out", str(tmp_path / "x.txt"))
    assert r.returncode == 2
    assert "feature_types_present" in r.stderr
    assert not (tmp_path / "x.txt").exists(), "aborted run still wrote output"


def test_low_id_overlap_aborts_with_exit_code_2(cli, tmp_path: Path, mini_gff: Path) -> None:
    bad = tmp_path / "XvsY_intersect.txt"
    bad.write_text("ID\tExpression\nNOT_IN_GFF_1\tUp_X_Down_Y\nNOT_IN_GFF_2\tUp_X_Down_Y\n")
    lnc = tmp_path / "lnc.txt"
    lnc.write_text("XR_001.1\n")
    r = cli("exclusive", "--lnc-ids", str(lnc), "--outdir", str(tmp_path),
            "--gff", str(mini_gff), "--level", "gene", "--intersects", str(bad))
    assert r.returncode == 2
    assert "expression_ids_in_gff" in r.stderr


def test_no_validate_skips_the_contract(cli, tmp_path: Path, fixtures: Path) -> None:
    r = cli("extract-ids", "--gff", str(fixtures / "ensembl_style.gff"),
            "--out", str(tmp_path / "x.txt"), "--no-validate")
    assert r.returncode == 0, r.stderr
    assert "validation skipped" in r.stderr


def test_drop_totals_are_printed_at_the_end(cli, tmp_path: Path, fixtures: Path) -> None:
    r = cli("extract-ids", "--gff", str(fixtures / "malformed.gff"), "--out", str(tmp_path / "x.txt"))
    assert "records dropped: 4" in r.stderr
    assert "non_integer_coordinates" in r.stderr
