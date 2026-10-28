"""Kernel helpers for the protein-inspector skill.

Everything of substance lives in the `protein_inspector` package, so a
notebook, a CLI, an MCP client and an agent all reach the same code. This file
is the Claude Science adapter over it: find the engine, then forward.

The wrappers take `*args, **kwargs` rather than restating each signature --
a copied signature here is a signature that drifts from the package's. The
real parameters and docstrings are on `protein_inspector.<name>`, and SKILL.md
documents the ones an agent needs.
"""
import os
import sys

SEARCH_PATHS = ("~/code/protein-inspector", "~/protein-inspector")

# THE OLDEST ENGINE THIS SKILL'S TEXT IS TRUE OF. SKILL.md documents layer
# descriptions and partners that morph with the target, which 1.13.0 does not
# have. Without a floor, `pip install protein-inspector` leaves an older release
# in place (the requirement is "already satisfied") and the agent renders with
# an engine that quietly ignores half of what it was told to pass.
MIN_VERSION = (1, 15, 0)


def inspector_engine(path=None):
    """Import the renderer, adding a local checkout to sys.path if needed."""
    candidates = []
    if path:
        candidates.append(path)
    if os.environ.get("PROTEIN_INSPECTOR_HOME"):
        candidates.append(os.environ["PROTEIN_INSPECTOR_HOME"])
    candidates.extend(os.path.expanduser(p) for p in SEARCH_PATHS)
    # FIRST MATCH WINS, and it has to actually contain the package. Inserting
    # every hit at sys.path[0] in turn reversed the order, so the hard-coded
    # fallback outranked an explicit PROTEIN_INSPECTOR_HOME and loaded the
    # wrong checkout without saying so.
    for candidate in candidates:
        if candidate and os.path.isdir(os.path.join(candidate, "protein_inspector")):
            if candidate not in sys.path:
                sys.path.insert(0, candidate)
            break
    try:
        import protein_inspector as engine
    except ImportError as exc:
        raise ImportError(
            "protein_inspector not importable. Either "
            "`pip install protein-inspector` or point PROTEIN_INSPECTOR_HOME "
            "at a local checkout."
        ) from exc
    have = engine_version(engine)
    if have < MIN_VERSION:
        raise ImportError(
            "protein_inspector %s is older than this skill needs (%s). Run "
            "`pip install -U \"protein-inspector>=%s\"` and restart the kernel "
            "(an already-imported package stays the old one until then), or "
            "point PROTEIN_INSPECTOR_HOME at a current checkout."
            % (".".join(map(str, have)), ".".join(map(str, MIN_VERSION)),
               ".".join(map(str, MIN_VERSION[:2])))
        )
    return engine


def engine_version(engine):
    """(major, minor, patch) of an imported engine, read off INSPECTOR_VERSION.

    That slug (`protein-inspector-1.15.0`) is what every bundle is stamped
    with, so it is the number a reader can match against what they installed.
    An engine with no readable version counts as (0, 0, 0): too old to trust.
    """
    import re
    slug = str(getattr(engine, "INSPECTOR_VERSION", ""))
    nums = [int(n) for n in re.findall(r"\d+", slug.rsplit("-", 1)[-1])][:3]
    return tuple(nums) + (0,) * (3 - len(nums))


def inspect_structure(*args, **kwargs):
    """Render one self-contained annotated HTML bundle. Returns the report.

    See `protein_inspector.inspect_structure` for the full signature.
    """
    return inspector_engine().inspect_structure(*args, **kwargs)


def inspection_table(*args, **kwargs):
    """Flat per-residue rows (chain, number, colour, layers, labels) from a bundle."""
    return inspector_engine().inspection_table(*args, **kwargs)


def build_manifest(*args, **kwargs):
    """Turn simple layer dicts into a validated inspection manifest."""
    return inspector_engine().build_manifest(*args, **kwargs)


def residue_order(*args, **kwargs):
    """[(chain_id, author_resnum)] in file order -- the canonical_position basis."""
    return inspector_engine().residue_order(*args, **kwargs)


def to_mmcif(*args, **kwargs):
    """Return a path to a .cif for `structure`, converting a .pdb if needed."""
    return inspector_engine().to_mmcif(*args, **kwargs)


def capture_bundle(*args, **kwargs):
    """Write a bundle's picture to a transparent PNG, for a report or a slide.

    Needs the [capture] extra and a browser. See
    `protein_inspector.capture.capture_bundle`.
    """
    return inspector_engine().capture_bundle(*args, **kwargs)


def default_palette():
    """The layer colours used when a layer does not name one."""
    return inspector_engine().DEFAULT_PALETTE
