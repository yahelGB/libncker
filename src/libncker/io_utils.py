from __future__ import annotations

from pathlib import Path
from typing import Iterator


def read_ids_one_per_line(path: Path) -> set[str]:
    ids: set[str] = set()
    with path.open("r", encoding="utf-8", errors="replace") as f:
        for line in f:
            s = line.strip()
            if not s:
                continue
            if s.lower() in {"id", "ids", "lncrna", "lncrna_ids"}:
                continue
            ids.add(s)
    return ids


def iter_tsv_rows(path: Path) -> Iterator[list[str]]:
    with path.open("r", encoding="utf-8", errors="replace") as f:
        for line in f:
            line = line.rstrip("\n")
            if not line:
                continue
            yield line.split("\t")


def detect_id_and_expression_columns(header: list[str]) -> tuple[int, int] | None:
    norm = [h.strip() for h in header]
    try:
        id_idx = norm.index("ID")
    except ValueError:
        return None
    try:
        expr_idx = norm.index("Expression")
    except ValueError:
        return None
    return id_idx, expr_idx


def top_prefixes(ids: set[str], n: int = 10) -> dict[str, int]:
    counts: dict[str, int] = {}
    for x in ids:
        if "_" in x:
            pref = x.split("_", 1)[0]
        elif "-" in x:
            pref = x.split("-", 1)[0]
        else:
            pref = x[:4]
        counts[pref] = counts.get(pref, 0) + 1
    return dict(sorted(counts.items(), key=lambda kv: -kv[1])[:n])