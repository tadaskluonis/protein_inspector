"""Render and read inspection bundles from the shell.

    protein-inspector render spec.json                 # spec -> bundle
    protein-inspector read inspection.html             # whole state as JSON
    protein-inspector read inspection.html residues    # just the residue table
    protein-inspector capture inspection.html fig.png  # the picture, as a PNG

`capture` needs a browser and is an optional extra:
`pip install "protein-inspector[capture]"` then `playwright install chromium`.
It writes a transparent PNG at 300 dpi by default; --dpi and --chrome override.

A render spec is JSON with the same keys `inspect_structure` takes:

    {"structure": "target.cif",
     "out": "inspection.html",
     "layers": [{"id": "hot", "label": "Hotspots", "color": "#dc2626",
                 "residues": {"56": "largest effect", "60": "second"}}],
     "about": {"What this is": "..."},
     "partner_chains": ["B"],
     "conformers": [{"path": "open.cif", "label": "open"}]}
"""
from __future__ import annotations

import json
import sys


def main(argv: list[str]) -> int:
    if not argv or argv[0] in {"-h", "--help"}:
        print(__doc__)
        return 0

    verb, rest = argv[0], argv[1:]

    if verb == "render":
        if not rest:
            print("render needs a spec.json", file=sys.stderr)
            return 2
        from .inspect import inspect_structure
        spec = json.loads(open(rest[0], encoding="utf-8").read())
        structure = spec.pop("structure")
        report = inspect_structure(structure, **spec)
        # The path first and alone on stdout, so a caller can pipe it.
        print(report["path"])
        if report.get("unmapped_residues"):
            print("unmapped residues (not in the model, so not painted): %r"
                  % (report["unmapped_residues"],), file=sys.stderr)
        return 0

    if verb == "capture":
        if len(rest) < 2:
            print("capture needs a bundle and an output .png", file=sys.stderr)
            return 2
        from .capture import capture_bundle
        bundle, out = rest[0], rest[1]
        opts = {}
        for flag, key, cast in (("--dpi", "dpi", int), ("--chrome", "chrome", str),
                                ("--settle-ms", "settle_ms", int)):
            if flag in rest:
                opts[key] = cast(rest[rest.index(flag) + 1])
        report = capture_bundle(bundle, out, **opts)
        print(report["path"])
        print("%(width)sx%(height)s, alpha=%(has_alpha)s, %(bytes)s bytes"
              % report, file=sys.stderr)
        return 0

    if verb == "read":
        from .renderer import read_inspection_bundle
        if not rest:
            print("read needs a bundle path", file=sys.stderr)
            return 2
        state = read_inspection_bundle(rest[0])
        if len(rest) > 1:
            if rest[1] not in state:
                print("no such key: %s; have %s" % (rest[1], sorted(state)),
                      file=sys.stderr)
                return 2
            state = state[rest[1]]
        json.dump(state, sys.stdout, indent=2)
        sys.stdout.write("\n")
        return 0

    # A bare path still reads, which is what the old entry point did.
    from .renderer import read_inspection_bundle
    state = read_inspection_bundle(verb)
    if rest:
        state = state.get(rest[0], state)
    json.dump(state, sys.stdout, indent=2)
    sys.stdout.write("\n")
    return 0


def cli() -> int:
    try:
        return main(sys.argv[1:])
    except BrokenPipeError:
        # `protein-inspector read x.html residues | head` closes the pipe under
        # us. Python would otherwise print a traceback on the way out, which
        # in a shell pipeline looks like the tool failed.
        try:
            sys.stdout.close()
        except BrokenPipeError:
            pass
        return 0


if __name__ == "__main__":
    raise SystemExit(cli())
