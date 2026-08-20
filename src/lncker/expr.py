from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True)
class ParsedExpression:
    up: str
    down: str


def parse_expression(expr: str) -> ParsedExpression | None:
    """
    Supports IDEAMEX patterns:
      - Down_<A>_Up_<B>
      - Up_<A>_Down_<B>

    Returns ParsedExpression(up=<tissue_up>, down=<tissue_down>) or None if invalid.
    Robust to underscores inside tissue names because we split by '_Up_' / '_Down_' tokens.
    """
    expr = expr.strip()
    if not expr:
        return None

    if expr.startswith("Down_"):
        if "_Up_" not in expr:
            return None
        left, right = expr.split("_Up_", 1)
        down = left[len("Down_") :].strip()
        up = right.strip()
        if not up or not down:
            return None
        return ParsedExpression(up=up, down=down)

    if expr.startswith("Up_"):
        if "_Down_" not in expr:
            return None
        left, right = expr.split("_Down_", 1)
        up = left[len("Up_") :].strip()
        down = right.strip()
        if not up or not down:
            return None
        return ParsedExpression(up=up, down=down)

    return None