"""Chain colouring: the palette, the per-chain overrides, and the config path.

WHY THIS FILE EXISTS. `chain_palette` was plumbed from `py2Dmol.view()` into
`config.color` and read back out in `parts/ui.js`, and the emitted JSON was
checked and found to carry it -- but `normalizeConfig` in core/mol.js rebuilds
`color` from a NAMED WHITELIST, so the key was dropped in between and every
chain still drew bright green. Asserting that a value appears in the emitted
page says nothing about whether the code that consumes it kept it.

So these tests EXECUTE the real src/core/mol.js rather than grepping it. The
system node is 12.x and cannot parse the optional chaining the source uses;
QuickJS (ES2020) can, and the functions under test are pure.
"""

import hashlib
import json
import pathlib
import re

import pytest

quickjs = pytest.importorskip("quickjs", reason="QuickJS runs the real mol.js")

REPO = pathlib.Path(__file__).resolve().parent.parent
MOL_JS = REPO / "src" / "core" / "mol.js"

# Enough of a browser for mol.js to finish evaluating its top level. Nothing
# under test touches the DOM; this only stops the module from throwing on load.
BROWSER_STUB = """
var window = globalThis; var self = globalThis;
var navigator = {userAgent: 'quickjs', gpu: undefined};
function _el() {
    return {style: {}, children: [],
            classList: {add() {}, remove() {}, contains() { return false; }},
            appendChild() {}, setAttribute() {}, addEventListener() {},
            querySelector() { return null; }, querySelectorAll() { return []; },
            getContext() { return null; }};
}
var document = {createElement: _el, createElementNS: _el, createDocumentFragment: _el,
                body: _el(), getElementById() { return null; },
                querySelector() { return null; }, querySelectorAll() { return []; },
                addEventListener() {}, dispatchEvent() {}};
var requestAnimationFrame = function () { return 0; };
var cancelAnimationFrame = function () {};
var CustomEvent = function (n, o) { this.type = n; Object.assign(this, o || {}); };
var performance = {now: function () { return 0; }};
"""


@pytest.fixture(scope="module")
def js():
    ctx = quickjs.Context()
    ctx.eval(BROWSER_STUB)
    ctx.eval(MOL_JS.read_text())
    return lambda expr: json.loads(ctx.eval("JSON.stringify(" + expr + ")"))


# --- the regression ------------------------------------------------------

def test_chain_palette_survives_config_normalization(js):
    """The whitelist in normalizeConfig must name chain_palette.

    This is the bug. Without the key in that object literal the palette is
    dropped silently between the page and the renderer -- no error, no warning,
    just the default palette, whose first colour is bright green.
    """
    out = js("normalizeConfig({color: {mode: 'chain', chain_palette: 'greys'}}).color")
    assert out["chain_palette"] == "greys"
    assert out["mode"] == "chain"


def test_chain_palette_accepts_the_legacy_flat_form(js):
    """ss_palette honours a top-level key as well as a nested one; so must this."""
    assert js("normalizeConfig({chain_palette: 'greys'}).color.chain_palette") == "greys"


def test_chain_palette_is_absent_rather_than_guessed_when_unset(js):
    """Unset must stay unset: ui.js only assigns when the key is truthy, and a
    default invented here would override a palette set directly on the renderer."""
    assert js("normalizeConfig({}).color.chain_palette === undefined") is True


# --- palette selection ---------------------------------------------------

def test_greys_palette_is_actually_grey(js):
    """Every entry is neutral (r == g == b) and neither pure white nor black."""
    greys = js("chainPaletteFor({chainPalette: 'greys'})")
    assert len(greys) >= 8, "too few greys and chains start repeating colours"
    for hex_color in greys:
        assert re.fullmatch(r"#[0-9a-f]{6}", hex_color), hex_color
        r, g, b = (int(hex_color[i:i + 2], 16) for i in (1, 3, 5))
        assert r == g == b, f"{hex_color} is not neutral"
        assert 0x30 <= r <= 0xe8, f"{hex_color} is too close to the background or to black"
    assert len(set(greys)) == len(greys), "duplicate greys make two chains indistinguishable"


def test_adjacent_chains_are_told_apart(js):
    """The ramp is interleaved light/dark on purpose: a monotonic ramp puts the
    two nearest values on the two chains most likely to be side by side."""
    greys = js("chainPaletteFor({chainPalette: 'greys'})")
    vals = [int(h[1:3], 16) for h in greys]
    gaps = [abs(b - a) for a, b in zip(vals, vals[1:])]
    assert min(gaps) >= 0x20, f"adjacent greys too close: {gaps}"


def test_default_palette_is_unchanged(js):
    """The greys are opt-in. Upstream py2Dmol must colour chains as it always
    did for anything that does not ask -- including an unknown palette name."""
    default = js("chainPaletteFor({})")
    assert default[0] == "#66ff66"
    assert js("chainPaletteFor({chainPalette: 'nope'})") == default
    assert js("chainPaletteFor({chainPalette: null})") == default


def test_colorblind_mode_outranks_the_palette(js):
    """Colourblind-safe is an accessibility setting, not a preference, so it
    wins over a palette choice rather than being silently overridden by one."""
    cb = js("chainPaletteFor({chainPalette: 'greys', colorblindMode: true})")
    assert cb != js("chainPaletteFor({chainPalette: 'greys'})")
    assert cb == js("chainPaletteFor({colorblindMode: true})")


def test_palette_names_are_discoverable(js):
    assert sorted(js("window.py2dmol_chainPalettes()")) == ["greys", "pymol"]


# --- per-chain overrides -------------------------------------------------

def test_override_beats_the_palette(js):
    r = "{chainPalette: 'greys', chainColorOverrides: {'A': '#ff0000'}}"
    assert js(f"chainColorHexFor({r}, 'A', 1)") == "#ff0000"


def test_override_applies_to_one_chain_only(js):
    """Keyed by chain, not global: pinning A must leave B on the palette."""
    r = "{chainPalette: 'greys', chainColorOverrides: {'A': '#ff0000'}}"
    assert js(f"chainColorHexFor({r}, 'B', 1)") == js("chainPaletteFor({chainPalette: 'greys'})[1]")


def test_override_is_revocable(js):
    """Deleting the key returns the chain to the palette -- the double-click
    reset in the panel does exactly this and nothing else."""
    palette_color = js("chainPaletteFor({chainPalette: 'greys'})[1]")
    assert js("chainColorHexFor({chainPalette: 'greys', chainColorOverrides: {}}, 'A', 1)") \
        == palette_color


def test_palette_wraps_rather_than_running_off_the_end(js):
    n = len(js("chainPaletteFor({chainPalette: 'greys'})"))
    assert js(f"chainColorHexFor({{chainPalette: 'greys'}}, 'Z', {n + 2})") \
        == js("chainColorHexFor({chainPalette: 'greys'}, 'Z', 2)")
