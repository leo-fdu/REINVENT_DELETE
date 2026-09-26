"""Pure-Python helpers to derive AutoDock Vina grid boxes from crystal ligands.

The docking box is defined once per target and shared by both generation
models so the comparison stays fair:

* center = heavy-atom centroid of the crystallographic ligand
* size   = ligand bounding box plus ``padding`` Angstrom on *each* side

This module has no third-party dependencies so it can be unit-tested
without RDKit/Meeko/Vina installed.
"""

from __future__ import annotations

from pathlib import Path

SECTION_PREFIX = "@<TRIPOS>"
ATOM_SECTION = "@<TRIPOS>ATOM"


def read_mol2_heavy_atom_coords(path: str | Path) -> list[tuple[float, float, float]]:
    """Return heavy-atom (non-hydrogen) coordinates from a Tripos .mol2 file.

    Parses only the ``@<TRIPOS>ATOM`` section. The element is taken from the
    Tripos atom type (the part before the dot, e.g. ``C.ar`` -> ``C``).
    Raises ValueError on malformed records or when no heavy atoms are found.
    """
    path = Path(path)
    coords: list[tuple[float, float, float]] = []
    in_atom_section = False
    with path.open() as handle:
        for line in handle:
            stripped = line.strip()
            if stripped.startswith(SECTION_PREFIX):
                in_atom_section = stripped.upper() == ATOM_SECTION
                continue
            if not in_atom_section or not stripped:
                continue
            fields = stripped.split()
            if len(fields) < 6:
                raise ValueError(f"malformed ATOM record in {path}: {line!r}")
            element = fields[5].split(".", 1)[0].upper()
            if element == "H":
                continue
            x, y, z = float(fields[2]), float(fields[3]), float(fields[4])
            coords.append((x, y, z))
    if not coords:
        raise ValueError(f"no heavy atoms found in ATOM section of {path}")
    return coords


def box_from_coords(
    coords: list[tuple[float, float, float]],
    padding: float = 8.0,
    ndigits: int = 3,
) -> dict:
    """Compute a Vina grid box from ligand coordinates.

    Returns ``{"center": [cx, cy, cz], "size": [sx, sy, sz]}`` where the
    center is the heavy-atom centroid and each size axis is the bounding-box
    extent plus ``padding`` on both sides (i.e. extent + 2 * padding).
    """
    xs, ys, zs = zip(*coords)
    n = len(coords)
    center = [round(sum(axis) / n, ndigits) for axis in (xs, ys, zs)]
    size = [round(max(axis) - min(axis) + 2.0 * padding, ndigits) for axis in (xs, ys, zs)]
    return {"center": center, "size": size}


def box_from_mol2(path: str | Path, padding: float = 8.0) -> dict:
    """Compute the Vina grid box for the crystal ligand stored in a .mol2 file."""
    return box_from_coords(read_mol2_heavy_atom_coords(path), padding=padding)
