"""Local-only structural inspection bundle renderer.

Two surfaces, same engine:

    inspect_structure(...)        layers in, one HTML bundle out -- start here
    render_inspection_bundle(...) the strict manifest schema underneath it

and `inspection_table(path)` / `read_inspection_bundle(path)` to read a
rendered bundle back and check that the annotation landed where you meant.
"""

from .renderer import (INSPECTOR_VERSION, read_inspection_bundle,
                       render_inspection_bundle)
from .inspect import (DEFAULT_PALETTE, build_manifest, inspect_structure,
                      inspection_table, residue_order, to_mmcif)

__all__ = [
    "INSPECTOR_VERSION",
    "inspect_structure",
    "inspection_table",
    "render_inspection_bundle",
    "read_inspection_bundle",
    "build_manifest",
    "residue_order",
    "to_mmcif",
    "DEFAULT_PALETTE",
]
