"""Expose the inspector to any MCP client.

    uvx --from "protein-inspector[mcp]" protein-inspector-mcp

Nothing here knows which agent is calling. The tool returns the PATH of a
self-contained HTML file; the host decides whether to open it, attach it or
hand it to a human.
"""
from __future__ import annotations

from typing import Any

try:
    from mcp.server.fastmcp import FastMCP
except ImportError as exc:  # pragma: no cover - import-time guidance only
    raise ImportError(
        "the MCP server needs the optional dependency: "
        "pip install 'protein-inspector[mcp]'"
    ) from exc

from .inspect import inspect_structure, inspection_table

mcp = FastMCP("protein-inspector")


@mcp.tool()
def render_bundle(
    structure: str,
    layers: list[dict[str, Any]],
    out: str = "inspection.html",
    chain: str | None = None,
    about: dict[str, str] | None = None,
    partner_chains: list[str] | None = None,
    conformers: list[dict[str, Any]] | None = None,
) -> dict[str, Any]:
    """Render one self-contained interactive HTML page for a structure.

    `layers` is a list of {"id", "label", "color", "residues"} where residues
    is a list of author residue numbers or a {number: note} mapping. Show few
    layers: a page with every layer you could compute is a legend with a
    structure behind it.

    Returns {"path", "unmapped_residues", "size_warning"}. Always read
    unmapped_residues -- those residues are not in the model and were not
    painted.
    """
    report = inspect_structure(
        structure, layers=layers, out=out, chain=chain, about=about,
        partner_chains=partner_chains, conformers=conformers)
    return {
        "path": report["path"],
        "unmapped_residues": report.get("unmapped_residues", []),
        "size_warning": report.get("size_warning"),
    }


@mcp.tool()
def read_bundle(path: str) -> list[dict[str, Any]]:
    """One row per modelled residue in a rendered bundle.

    Chain, author number, residue name, coordinates, colour, layers, labels --
    what the reader will actually see, which is how you check that an
    annotation landed on the residue you named.
    """
    return inspection_table(path)


def main() -> None:
    mcp.run()


if __name__ == "__main__":
    main()
