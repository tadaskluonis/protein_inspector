"""Dump a rendered bundle's inspection state as JSON, for callers that would
rather not import anything.

    python -m bindos_structure_inspector inspection.html            # whole state
    python -m bindos_structure_inspector inspection.html residues   # just the table
"""
from __future__ import annotations

import json
import sys

from .renderer import read_inspection_bundle


def main(argv: list[str]) -> int:
    if not argv or argv[0] in {"-h", "--help"}:
        print(__doc__)
        return 0
    state = read_inspection_bundle(argv[0])
    if len(argv) > 1:
        if argv[1] not in state:
            print(f"no such key: {argv[1]}; have {sorted(state)}", file=sys.stderr)
            return 2
        state = state[argv[1]]
    json.dump(state, sys.stdout, indent=2)
    sys.stdout.write("\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
