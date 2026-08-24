"""
Provenance headers for every lncker output.

Each output file opens with a block of ``#`` comment lines recording exactly
what produced it: tool version, git SHA, md5 of every input, the full command
line, the RNG seed, and a UTC timestamp. Readers in :mod:`lncker.io_utils`
skip ``#`` lines, so the header is invisible to the rest of the pipeline.

Note for downstream consumers (R, pandas): read these files with a comment
character set, e.g. ``read.delim(f, comment.char = "#")`` or
``pd.read_csv(f, sep="\\t", comment="#")``.
"""

from __future__ import annotations

import hashlib
import os
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence

from lncker import __version__

PROVENANCE_PREFIX = "#"

# lncker performs no random sampling; every command is a pure function of its
# inputs. The seed is recorded anyway so the header format stays stable if a
# future command does randomise.
DEFAULT_RNG_SEED = 0

_MD5_CHUNK = 1024 * 1024


def md5sum(path: Path) -> str:
    h = hashlib.md5()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(_MD5_CHUNK), b""):
            h.update(chunk)
    return h.hexdigest()


def git_sha(repo_hint: Optional[Path] = None) -> str:
    """Short git SHA of the working tree lncker was run from, or 'unknown'."""
    start = Path(repo_hint) if repo_hint else Path(__file__).resolve().parent
    try:
        res = subprocess.run(
            ["git", "-C", str(start), "rev-parse", "--short", "HEAD"],
            capture_output=True,
            text=True,
            timeout=10,
        )
    except (OSError, subprocess.SubprocessError):
        return "unknown"
    if res.returncode != 0:
        return "unknown"
    sha = res.stdout.strip()
    if not sha:
        return "unknown"
    try:
        dirty = subprocess.run(
            ["git", "-C", str(start), "status", "--porcelain"],
            capture_output=True,
            text=True,
            timeout=10,
        )
        if dirty.returncode == 0 and dirty.stdout.strip():
            sha += "-dirty"
    except (OSError, subprocess.SubprocessError):
        pass
    return sha


class Provenance:
    """Collects provenance for one lncker command and renders it as a header."""

    def __init__(
        self,
        command: str,
        inputs: Optional[Sequence[Path]] = None,
        params: Optional[Dict[str, object]] = None,
        rng_seed: int = DEFAULT_RNG_SEED,
        argv: Optional[Sequence[str]] = None,
    ) -> None:
        self.command = command
        self.params = dict(params or {})
        self.rng_seed = rng_seed
        self.argv = list(argv if argv is not None else sys.argv)
        self.timestamp = datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")
        self.tool_version = __version__
        self.git_sha = git_sha()
        self.input_md5s: List[str] = []
        for p in inputs or []:
            self.add_input(Path(p))

    def add_input(self, path: Path) -> None:
        """Record ``path``'s md5. Unreadable inputs are recorded as such, not skipped."""
        try:
            digest = md5sum(path)
            size = path.stat().st_size
            self.input_md5s.append(f"{digest}  {size}  {path}")
        except OSError as exc:
            self.input_md5s.append(f"UNREADABLE  -  {path} ({exc.strerror})")

    def header_lines(self) -> List[str]:
        p = PROVENANCE_PREFIX
        lines = [
            f"{p} lncker provenance",
            f"{p} tool_version: {self.tool_version}",
            f"{p} git_sha: {self.git_sha}",
            f"{p} command: {self.command}",
            f"{p} timestamp_utc: {self.timestamp}",
            f"{p} rng_seed: {self.rng_seed} (lncker is deterministic; no RNG is used)",
            f"{p} python: {sys.version.split()[0]}",
            f"{p} cwd: {os.getcwd()}",
            f"{p} argv: {' '.join(self.argv)}",
        ]
        if self.params:
            rendered = " ".join(f"{k}={v}" for k, v in sorted(self.params.items()))
            lines.append(f"{p} params: {rendered}")
        if self.input_md5s:
            lines.append(f"{p} inputs (md5  bytes  path):")
            lines.extend(f"{p}   {row}" for row in self.input_md5s)
        return lines

    def write_header(self, fh) -> None:
        for line in self.header_lines():
            fh.write(line + "\n")


def write_gff_provenance_file(
    gff: Path,
    out: Optional[Path] = None,
    accession: str = "",
    annotation_release: str = "",
    download_date: str = "",
    extra: Optional[Dict[str, str]] = None,
) -> Path:
    """Write ``PROVENANCE.txt`` next to the reference GFF.

    Accession and annotation release are taken from the GFF's own NCBI pragmas
    when the caller does not supply them.
    """
    from lncker.validate import read_gff_provenance_pragmas

    pragmas = read_gff_provenance_pragmas(gff)
    accession = accession or pragmas.get("genome-build-accession", "")
    annotation_release = annotation_release or pragmas.get("annotation-source", "")
    target = out or gff.parent / "PROVENANCE.txt"

    rows: List[str] = [
        "# Reference annotation provenance",
        f"file: {gff.name}",
        f"bytes: {gff.stat().st_size}",
        f"md5: {md5sum(gff)}",
        f"accession: {accession or 'UNKNOWN'}",
        f"annotation_release: {annotation_release or 'UNKNOWN'}",
        f"genome_build: {pragmas.get('genome-build', 'UNKNOWN')}",
        f"gff_version: {pragmas.get('gff-version', 'UNKNOWN')}",
        f"processor: {pragmas.get('processor', 'UNKNOWN')}",
        f"download_date: {download_date or 'UNKNOWN'}",
        f"recorded_utc: {datetime.now(timezone.utc).strftime('%Y-%m-%dT%H:%M:%SZ')}",
        f"recorded_by: lncker {__version__} ({git_sha()})",
    ]
    for k, v in sorted((extra or {}).items()):
        rows.append(f"{k}: {v}")

    target.write_text("\n".join(rows) + "\n", encoding="utf-8")
    return target
