"""Conformational morphs: residue mapping, superposition, interpolation.

The morph is a DEPICTION OF ITS ENDPOINTS, not a pathway -- intermediates are
Cartesian interpolations, so they are not physical and bond geometry is not
preserved. Everything here tests that the endpoints are exactly the structures
they claim to be and that the residues paired between them are the right ones.

The residue-mapping guard carries most of the weight. A morph built over a
quiet intersection of two residue sets looks exactly as smooth as a correct
one while interpolating the wrong pairs, so a mismatch has to be refused
loudly rather than worked around.
"""

import hashlib
import json
import pathlib

import numpy as np
import pytest
from Bio.PDB import Atom, Chain, Model, Residue, Structure
from Bio.PDB.mmcifio import MMCIFIO

from protein_inspector import render_inspection_bundle
from protein_inspector.renderer import (_ca_trace, _conformer_stacks, _morph_frames,
                                        _superpose)


def _write(path, coords, chain_id="A", numbers=None, resname="ALA", extra_atom=False):
    """A single-chain CA-only structure at the given coordinates."""
    numbers = numbers if numbers is not None else range(1, len(coords) + 1)
    structure = Structure.Structure("s")
    model = Model.Model(0)
    chain = Chain.Chain(chain_id)
    for number, xyz in zip(numbers, coords):
        residue = Residue.Residue((" ", int(number), " "), resname, " ")
        residue.add(Atom.Atom("CA", tuple(float(v) for v in xyz), 0.0, 1.0, " ",
                              "CA", int(number), element="C"))
        if extra_atom:
            residue.add(Atom.Atom("CB", tuple(float(v) + 1.0 for v in xyz), 0.0, 1.0,
                                  " ", "CB", int(number) + 1000, element="C"))
        chain.add(residue)
    model.add(chain)
    structure.add(model)
    io = MMCIFIO()
    io.set_structure(structure)
    io.save(str(path))
    return path


def _helix(n=12, twist=0.6, rise=1.5, radius=4.0, phase=0.0):
    t = np.arange(n) * twist + phase
    return np.stack([radius * np.cos(t), radius * np.sin(t), np.arange(n) * rise], 1)


# --- residue mapping -----------------------------------------------------

def test_ca_trace_keys_on_chain_and_author_number(tmp_path):
    trace = _ca_trace(_write(tmp_path / "a.cif", _helix(4), chain_id="B",
                             numbers=[10, 11, 12, 13]))
    assert sorted(trace) == [("B", 10), ("B", 11), ("B", 12), ("B", 13)]
    assert all(name == "ALA" and len(xyz) == 3 for name, xyz in trace.values())


def test_ca_trace_takes_one_atom_per_residue(tmp_path):
    """Sidechain atoms must not enter the trace, or the frames carry more
    points than the colour array has residues."""
    assert len(_ca_trace(_write(tmp_path / "a.cif", _helix(5), extra_atom=True))) == 5


def test_ca_trace_refuses_a_duplicate_residue_address(tmp_path):
    """An insertion code or altloc collapsing two residues onto one address
    makes the pairing ambiguous, and ambiguity here is silent."""
    path = tmp_path / "dup.cif"
    coords = _helix(3)
    structure = Structure.Structure("s")
    model = Model.Model(0)
    chain = Chain.Chain("A")
    # 3 and 3A are different residues to Biopython and the same address to us.
    for number, icode, xyz in [(1, " ", coords[0]), (3, " ", coords[1]),
                               (3, "A", coords[2])]:
        residue = Residue.Residue((" ", number, icode), "ALA", " ")
        residue.add(Atom.Atom("CA", tuple(float(v) for v in xyz), 0.0, 1.0, " ",
                              "CA", 1, element="C"))
        chain.add(residue)
    model.add(chain)
    structure.add(model)
    io = MMCIFIO()
    io.set_structure(structure)
    io.save(str(path))
    with pytest.raises(ValueError, match="twice|ambiguous"):
        _ca_trace(path)


def test_morph_refuses_a_mismatched_residue_set_and_names_the_offenders(tmp_path):
    ref = _write(tmp_path / "ref.cif", _helix(8), numbers=range(1, 9))
    short = _write(tmp_path / "short.cif", _helix(6, phase=0.3), numbers=range(1, 7))
    with pytest.raises(ValueError) as excinfo:
        _morph_frames(ref, [{"path": str(short)}], steps=4)
    message = str(excinfo.value)
    assert "not exact" in message
    # The specific residues, not just a count: "they differ" is not actionable.
    assert "('A', 7)" in message and "('A', 8)" in message


def test_morph_refuses_a_conformer_whose_hash_does_not_match(tmp_path):
    ref = _write(tmp_path / "ref.cif", _helix(6))
    other = _write(tmp_path / "other.cif", _helix(6, phase=0.4))
    with pytest.raises(ValueError, match="hash does not match"):
        _morph_frames(ref, [{"path": str(other), "sha256": "0" * 64}], steps=2)


def test_morph_accepts_a_matching_hash(tmp_path):
    ref = _write(tmp_path / "ref.cif", _helix(6))
    other = _write(tmp_path / "other.cif", _helix(6, phase=0.4))
    digest = hashlib.sha256(other.read_bytes()).hexdigest()
    _, _, frames, _ = _morph_frames(ref, [{"path": str(other), "sha256": digest}], steps=2)
    assert len(frames) == 3


# --- superposition -------------------------------------------------------

def test_superpose_recovers_a_rigid_transform(tmp_path):
    target = _helix(20)
    angle = 0.7
    rotation = np.array([[np.cos(angle), -np.sin(angle), 0],
                         [np.sin(angle), np.cos(angle), 0], [0, 0, 1.0]])
    moved = target @ rotation.T + np.array([13.0, -4.0, 7.5])
    assert np.allclose(_superpose(moved, target), target, atol=1e-6)


def test_superpose_forbids_reflection(tmp_path):
    """A mirrored chain is not the same chain. Kabsch without the determinant
    correction will happily flip one onto the other and report RMSD 0."""
    target = _helix(20)
    mirrored = target * np.array([1.0, 1.0, -1.0])
    rmsd = np.sqrt(((_superpose(mirrored, target) - target) ** 2).sum(1).mean())
    assert rmsd > 1.0, "reflection was allowed, so a mirror image superposes perfectly"


# --- interpolation -------------------------------------------------------

def test_endpoints_are_the_structures_not_interpolations(tmp_path):
    """First and last frame must be the inputs exactly (the last after
    superposition), or the morph never actually shows either conformer."""
    a, b = _helix(10), _helix(10, phase=0.9)
    ref = _write(tmp_path / "a.cif", a)
    other = _write(tmp_path / "b.cif", b)
    keys, names, frames, report = _morph_frames(ref, [{"path": str(other)}], steps=5)
    # mmCIF carries three decimals, so the round trip is exact to ~1e-3 A,
    # not to machine precision. The point is that the endpoint is the
    # structure rather than an interpolation of it.
    assert np.allclose(frames[0], a, atol=5e-3)
    assert np.allclose(frames[-1], _superpose(b, a), atol=5e-3)
    assert report["n_frames"] == len(frames) == 5 + 1
    assert report["n_residues"] == len(keys) == len(names) == 10


def test_interpolation_is_linear_in_time(tmp_path):
    a, b = _helix(8), _helix(8, phase=1.1)
    ref = _write(tmp_path / "a.cif", a)
    other = _write(tmp_path / "b.cif", b)
    _, _, frames, _ = _morph_frames(ref, [{"path": str(other)}], steps=2)
    assert np.allclose(frames[1], (frames[0] + frames[2]) / 2.0, atol=1e-6)


def test_multiple_conformers_chain_end_to_end(tmp_path):
    a, b, c = _helix(9), _helix(9, phase=0.5), _helix(9, phase=1.0)
    ref = _write(tmp_path / "a.cif", a)
    paths = [{"path": str(_write(tmp_path / "b.cif", b)), "label": "mid"},
             {"path": str(_write(tmp_path / "c.cif", c))}]
    _, _, frames, report = _morph_frames(ref, paths, steps=4)
    assert report["conformers"] == ["reference", "mid", "c"]
    assert report["n_frames"] == len(frames) == 4 * 2 + 1
    assert len(report["endpoint_rmsd_A"]) == 2


def test_reported_rmsd_is_the_superposed_endpoint_distance(tmp_path):
    a, b = _helix(15), _helix(15, phase=0.8)
    ref = _write(tmp_path / "a.cif", a)
    other = _write(tmp_path / "b.cif", b)
    _, _, _, report = _morph_frames(ref, [{"path": str(other)}], steps=3)
    expected = np.sqrt(((_superpose(b, a) - a) ** 2).sum(1).mean())
    assert report["endpoint_rmsd_A"] == [round(float(expected), 2)]


# --- the property the viewer depends on ----------------------------------

def test_morph_does_not_disturb_the_colouring(tmp_path):
    """THE POINT OF THE FEATURE. Colour is per-residue and frames carry only
    coordinates, so the annotation layers must paint identically with a morph
    attached and without one. If this fails, the morph is repainting."""
    a, b = _helix(6), _helix(6, phase=0.7)
    ref = _write(tmp_path / "a.cif", a)
    other = _write(tmp_path / "b.cif", b)
    manifest = {
        "schema_version": "protein-inspector-manifest-1",
        "annotations": [
            {"annotation_id": "p", "kind": "custom", "label": "probe",
             "color": "#dc2626",
             "residue": {"component_id": "target", "canonical_position": 3,
                         "chain_id": "A", "author_residue_number": 3},
             "resolved": True, "evidence_ids": ["e"], "method": "manual"},
        ],
    }
    digest = hashlib.sha256(ref.read_bytes()).hexdigest()

    def build(name, morph):
        render_inspection_bundle(
            mmcif_path=str(ref), mmcif_sha256=digest, inspection_manifest=manifest,
            output_dir=str(tmp_path / name), output_name="b", extras=True,
            **({"morph_to": morph, "morph_steps": 3} if morph else {}))
        state = json.loads((tmp_path / name / "b.viewer.json").read_text())
        return state["viewer"]["objects"][0]

    plain = build("plain", None)
    morphed = build("morphed", [{"path": str(other), "label": "open"}])

    assert morphed["color"] == plain["color"], "morph changed the colour assignment"
    assert len(plain["frames"]) == 1 and len(morphed["frames"]) == 4
    # ...and the extra frames are coordinates only: same residues, same order.
    assert morphed["frames"][0]["chains"] == plain["frames"][0]["chains"]
    assert (morphed["frames"][0]["residue_numbers"]
            == plain["frames"][0]["residue_numbers"])


# --- union mapping: states that differ by an indel -----------------------

def test_union_keeps_a_position_only_one_state_models(tmp_path):
    """THE POINT OF THE MODE. Under "intersection" a residue only one state has
    is removed from EVERY state, so the state whose file is complete draws a
    chain break at a position its own structure models -- which reads as a
    broken model rather than an alignment gap. "union" keeps the position."""
    full = _write(tmp_path / "full.cif", _helix(6), numbers=[1, 2, 3, 4, 5, 6])
    short = _write(tmp_path / "short.cif", _helix(5, phase=0.3),
                   numbers=[1, 2, 3, 5, 6])

    keys, names, stacks, labels, report, _, state_names = _conformer_stacks(
        full, [{"path": str(short), "label": "short"}], "union", "full")

    assert [k[1] for k in keys] == [1, 2, 3, 4, 5, 6]
    assert report["n_residues"] == 6 and report["n_shared_residues"] == 5
    assert report["conformers"][0]["residues_bridged"] == []
    assert report["conformers"][1]["residues_bridged"] == [["A", 4]]
    # ...and "intersection" is unchanged: the position goes, from both states.
    inter = _conformer_stacks(full, [{"path": str(short), "label": "short"}],
                              "intersection", "full")[4]
    assert inter["n_residues"] == 5


def test_union_names_the_missing_position_unk_in_the_state_that_lacks_it(tmp_path):
    """A BRIDGE MUST NOT CLAIM A RESIDUE. The strip and the residue card read
    the per-state names, so a state that does not model position 4 has to say
    so there rather than inherit the reference's residue."""
    full = _write(tmp_path / "full.cif", _helix(6), numbers=[1, 2, 3, 4, 5, 6])
    short = _write(tmp_path / "short.cif", _helix(5, phase=0.3),
                   numbers=[1, 2, 3, 5, 6])

    keys, _names, _stacks, _labels, _report, _fits, state_names = _conformer_stacks(
        full, [{"path": str(short), "label": "short"}], "union", "full")

    gap = [k[1] for k in keys].index(4)
    assert state_names[0][gap] == "ALA"
    assert state_names[1][gap] == "UNK"
    assert "UNK" not in state_names[0]


def test_union_bridges_on_the_lacking_states_own_path(tmp_path):
    """The bridge is a point between the flanks THAT STATE models, so its
    ribbon stays continuous along the bond that is really there -- it is not
    the reference's coordinate borrowed into another structure."""
    full = _write(tmp_path / "full.cif", _helix(6), numbers=[1, 2, 3, 4, 5, 6])
    short = _write(tmp_path / "short.cif", _helix(5, phase=0.3),
                   numbers=[1, 2, 3, 5, 6])

    keys, _n, stacks, _l, _r, _f, _sn = _conformer_stacks(
        full, [{"path": str(short), "label": "short"}], "union", "full")

    gap = [k[1] for k in keys].index(4)
    bridged = stacks[1][gap]
    midpoint = (stacks[1][gap - 1] + stacks[1][gap + 1]) / 2.0
    assert np.allclose(bridged, midpoint, atol=1e-9)


def test_union_keeps_bridges_out_of_the_fit_and_the_rmsd(tmp_path):
    """A DEPICTION MUST NOT MOVE THE SUPERPOSITION. Both numbers have to be the
    ones the shared positions alone would give, or an invented coordinate is
    quietly steering the comparison the reader is shown."""
    coords = _helix(6)
    full = _write(tmp_path / "full.cif", coords, numbers=[1, 2, 3, 4, 5, 6])
    moved = _helix(5, phase=0.3)
    short = _write(tmp_path / "short.cif", moved, numbers=[1, 2, 3, 5, 6])
    # the same five residues on both sides, so no bridge is needed at all
    five = _write(tmp_path / "five.cif", coords[[0, 1, 2, 4, 5]],
                  numbers=[1, 2, 3, 5, 6])

    union = _conformer_stacks(full, [{"path": str(short), "label": "s"}],
                              "union", "full")[4]
    exact = _conformer_stacks(five, [{"path": str(short), "label": "s"}],
                              "exact", "five")[4]
    assert (union["conformers"][1]["rmsd_to_reference_A"]
            == exact["conformers"][1]["rmsd_to_reference_A"])


def test_union_leaves_the_reference_frame_its_own_names(tmp_path):
    """REGRESSION. state_names[0] was the same list object the partner loop
    appends to, so the reference frame's name list grew by the partner tail
    twice and py2Dmol dropped names and types for that frame -- a warning on
    stdout and a frame that cannot label its own residues."""
    full = _write(tmp_path / "full.cif", _helix(6), numbers=[1, 2, 3, 4, 5, 6])
    short = _write(tmp_path / "short.cif", _helix(5, phase=0.3),
                   numbers=[1, 2, 3, 5, 6])

    keys, names, _s, _l, _r, _f, state_names = _conformer_stacks(
        full, [{"path": str(short), "label": "short"}], "union", "full")

    names.append("SENTINEL")
    assert len(state_names[0]) == len(keys)
    assert "SENTINEL" not in state_names[0]


def test_union_is_the_only_new_mode_and_a_bad_one_is_refused(tmp_path):
    full = _write(tmp_path / "full.cif", _helix(6))
    short = _write(tmp_path / "short.cif", _helix(6, phase=0.3))
    with pytest.raises(ValueError, match="union"):
        _conformer_stacks(full, [{"path": str(short), "label": "s"}],
                          "nearest", "full")


def test_an_altloc_on_a_backbone_atom_does_not_kill_the_trace(tmp_path):
    """Biopython hands back an Entity, and `a or b` asks it for a truth value.

    `Entity.__len__` is `len(self.selected_child)`, and on a DisorderedAtom
    whose selected child is a plain Atom that raises `TypeError: object of
    type 'Atom' has no len()`. Every crystal structure with an altloc at a CA
    therefore died on the way into a morph, with a traceback that named
    neither the file nor the residue. Found on PDB 2Y0Q / 2Y20 / 2Y21, the
    Af1503 A291 variant series, while building a five-state gearbox morph.
    """
    from pathlib import Path

    from Bio.PDB import Atom, Chain, Model, Residue, Structure
    from Bio.PDB.mmcifio import MMCIFIO

    from protein_inspector.renderer import _ca_trace

    structure = Structure.Structure("fixture")
    model = Model.Model(0)
    chain = Chain.Chain("A")
    for index in range(1, 6):
        residue = Residue.Residue((" ", index, " "), "ALA", " ")
        atom = Atom.Atom("CA", (float(index * 3), 0.0, 0.0), 0.0, 1.0, " ", "CA",
                         index, element="C")
        if index == 3:
            # A disordered CA carrying one altloc -- the shape that threw.
            atom.set_altloc("A")
            disordered = Atom.DisorderedAtom("CA")
            disordered.disordered_add(atom)
            residue.add(disordered)
        else:
            residue.add(atom)
        chain.add(residue)
    model.add(chain)
    structure.add(model)
    io = MMCIFIO()
    io.set_structure(structure)
    path = tmp_path / "altloc.cif"
    io.save(str(path))

    trace = _ca_trace(Path(str(path)))
    assert len(trace) == 5
    assert trace[("A", 3)][0] == "ALA"
