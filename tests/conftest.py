from __future__ import annotations

import subprocess
import sys
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parent.parent
SRC = REPO_ROOT / "src"
FIXTURES = Path(__file__).parent / "fixtures"

# Import lncker from the working tree, not from whatever is pip-installed.
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))


@pytest.fixture(scope="session")
def fixtures() -> Path:
    return FIXTURES


@pytest.fixture(scope="session")
def mini_gff(fixtures: Path) -> Path:
    return fixtures / "mini.gff"


@pytest.fixture(scope="session")
def intersects(fixtures: Path) -> list[Path]:
    return sorted(fixtures.glob("*vs*_intersect.txt"))


@pytest.fixture(scope="session")
def cli():
    """Run the lncker CLI as a subprocess and return the CompletedProcess."""

    def _run(*args: str, cwd: Path | None = None) -> subprocess.CompletedProcess:
        return subprocess.run(
            [sys.executable, "-m", "lncker", *args],
            capture_output=True,
            text=True,
            cwd=str(cwd or REPO_ROOT),
            env={"PYTHONPATH": str(SRC), "PATH": "/usr/bin:/bin", "HOME": str(Path.home())},
        )

    return _run
