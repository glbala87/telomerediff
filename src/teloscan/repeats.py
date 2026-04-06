"""Utilities for parsing k-mer specifications."""

from __future__ import annotations

from typing import List


def parse_k(k) -> List[int]:
    """
    Parse k-mer size specification into a sorted list of integers.

    Supported formats:
      - int or str of single int:  6   -> [6]
      - range:  "4-15"             -> [4, 5, ..., 15]
      - comma-separated:  "5,6,7"  -> [5, 6, 7]
      - mixed:  "4-6,8,10-12"     -> [4, 5, 6, 8, 10, 11, 12]

    Raises ValueError on invalid input.
    """
    if isinstance(k, int):
        if k < 1:
            raise ValueError(f"k must be >= 1, got {k}")
        return [k]

    if isinstance(k, (list, tuple)):
        result = sorted(set(int(x) for x in k))
        if any(v < 1 for v in result):
            raise ValueError(f"All k values must be >= 1, got {k}")
        return result

    # String parsing
    k_str = str(k).strip()
    if not k_str:
        raise ValueError("Empty k specification")

    values = set()
    for part in k_str.split(","):
        part = part.strip()
        if not part:
            continue
        if "-" in part:
            pieces = part.split("-")
            if len(pieces) != 2:
                raise ValueError(f"Invalid range: {part!r}")
            try:
                lo, hi = int(pieces[0].strip()), int(pieces[1].strip())
            except ValueError:
                raise ValueError(f"Non-integer in range: {part!r}")
            if lo > hi:
                raise ValueError(f"Invalid range {lo}-{hi}: start > end")
            values.update(range(lo, hi + 1))
        else:
            try:
                values.add(int(part))
            except ValueError:
                raise ValueError(f"Non-integer k value: {part!r}")

    result = sorted(values)
    if not result:
        raise ValueError(f"No valid k values parsed from: {k!r}")
    if any(v < 1 for v in result):
        raise ValueError(f"All k values must be >= 1, got values: {result}")
    return result
