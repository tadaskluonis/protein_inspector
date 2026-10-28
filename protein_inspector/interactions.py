"""Hydrogen bonds from a structure, as endpoint coordinates a viewer can draw.

Heavy-atom geometry only: structures from crystallography arrive without
hydrogens, so the criterion is the usual donor-acceptor distance plus a
direction test at the donor, taken against the heavy atom the donor is
covalently bonded to. No hydrogens are built and none are needed.
"""
from __future__ import annotations
import math
import gemmi

SC_DONORS = {
    "ARG": ("NE", "NH1", "NH2"), "ASN": ("ND2",), "GLN": ("NE2",),
    "HIS": ("ND1", "NE2"), "LYS": ("NZ",), "SER": ("OG",), "THR": ("OG1",),
    "TRP": ("NE1",), "TYR": ("OH",), "CYS": ("SG",),
}
SC_ACCEPTORS = {
    "ASP": ("OD1", "OD2"), "GLU": ("OE1", "OE2"), "ASN": ("OD1",),
    "GLN": ("OE1",), "HIS": ("ND1", "NE2"), "SER": ("OG",), "THR": ("OG1",),
    "TYR": ("OH",), "MET": ("SD",),
}
BB_DONOR = "N"
BB_ACCEPTORS = ("O", "OXT")


def _is_donor(res, atom):
    if atom.name == BB_DONOR:
        return res.name != "PRO"
    return atom.name in SC_DONORS.get(res.name, ())


def _is_acceptor(res, atom):
    if atom.name in BB_ACCEPTORS:
        return True
    return atom.name in SC_ACCEPTORS.get(res.name, ())


def _angle(a, b, c):
    """Angle a-b-c in degrees."""
    ux, uy, uz = a.x - b.x, a.y - b.y, a.z - b.z
    vx, vy, vz = c.x - b.x, c.y - b.y, c.z - b.z
    nu = math.sqrt(ux * ux + uy * uy + uz * uz)
    nv = math.sqrt(vx * vx + vy * vy + vz * vz)
    if nu == 0 or nv == 0:
        return 180.0
    cos = max(-1.0, min(1.0, (ux * vx + uy * vy + uz * vz) / (nu * nv)))
    return math.degrees(math.acos(cos))


def _antecedents(res, atom, cutoff=1.9):
    """The heavy atoms this one is bonded to, by distance, within its residue."""
    out = []
    for other in res:
        if other.name == atom.name or other.element == gemmi.Element("H"):
            continue
        if atom.pos.dist(other.pos) <= cutoff:
            out.append(other)
    return out


def find_hbonds(path, *, max_dist=3.5, min_angle=100.0, include_water=False,
                chains=None, model=0, skip_adjacent_backbone=True):
    st = gemmi.read_structure(str(path))
    st.remove_alternative_conformations()
    st.remove_hydrogens()
    if not include_water:
        st.remove_waters()
    st.setup_entities()
    mdl = st[model]

    sites = []                      # (chain_name, residue, atom)
    for ch in mdl:
        if chains and ch.name not in chains:
            continue
        for res in ch:
            for atom in res:
                sites.append((ch.name, res, atom))

    ns = gemmi.NeighborSearch(mdl, st.cell, max_dist + 1.0).populate()
    seen, bonds = set(), []
    for ch_name, res, atom in sites:
        if not _is_donor(res, atom):
            continue
        ante = _antecedents(res, atom)
        marks = ns.find_atoms(atom.pos, '\0', radius=max_dist)
        for mark in marks:
            cra = mark.to_cra(mdl)
            if chains and cra.chain.name not in chains:
                continue
            if cra.residue.seqid.num == res.seqid.num and cra.chain.name == ch_name:
                continue
            if not _is_acceptor(cra.residue, cra.atom):
                continue
            d = atom.pos.dist(cra.atom.pos)
            if d > max_dist:
                continue
            bb_pair = atom.name == BB_DONOR and cra.atom.name in BB_ACCEPTORS
            if (skip_adjacent_backbone and bb_pair and ch_name == cra.chain.name
                    and abs(cra.residue.seqid.num - res.seqid.num) < 2):
                continue
            if ante and min(_angle(a.pos, atom.pos, cra.atom.pos) for a in ante) < min_angle:
                continue
            key = (ch_name, res.seqid.num, atom.name,
                   cra.chain.name, cra.residue.seqid.num, cra.atom.name)
            if key in seen or key[3:] + key[:3] in seen:
                continue
            seen.add(key)
            d_bb = atom.name == BB_DONOR
            a_bb = cra.atom.name in BB_ACCEPTORS
            kind = "bb" if (d_bb and a_bb) else ("sc" if (not d_bb and not a_bb) else "mixed")
            bonds.append({
                "donor": {"chain": ch_name, "resnum": res.seqid.num,
                          "resname": res.name, "atom": atom.name,
                          "xyz": [round(atom.pos.x, 3), round(atom.pos.y, 3), round(atom.pos.z, 3)]},
                "acceptor": {"chain": cra.chain.name, "resnum": cra.residue.seqid.num,
                             "resname": cra.residue.name, "atom": cra.atom.name,
                             "xyz": [round(cra.atom.pos.x, 3), round(cra.atom.pos.y, 3),
                                     round(cra.atom.pos.z, 3)]},
                "distance": round(d, 2), "kind": kind,
            })
    bonds.sort(key=lambda b: (b["donor"]["chain"], b["donor"]["resnum"], b["acceptor"]["resnum"]))
    return bonds


def as_payload(bonds, colors=None):
    """The compact form the page plugin reads."""
    colors = colors or {"bb": "#3d7ebf", "mixed": "#3d7ebf", "sc": "#3d7ebf"}
    return {
        "bonds": [{
            "a": b["donor"]["xyz"], "b": b["acceptor"]["xyz"],
            "ac": b["donor"]["chain"], "ar": b["donor"]["resnum"],
            "bc": b["acceptor"]["chain"], "br": b["acceptor"]["resnum"],
            "k": b["kind"], "d": b["distance"],
        } for b in bonds],
        "colors": colors,
    }
