"""Shared cube file utilities: grid setup, file writing, and string generation."""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from chem_constants import ANG2BOHR, GRID_POINTS, GRID_MARGIN, ATOMIC_NUMBERS


def make_grid(atoms):
    """Create a grid around atom positions. atoms = [(sym,x,y,z), ...]."""
    coords = np.array([[x, y, z] for _, x, y, z in atoms]) * ANG2BOHR
    lo = coords.min(0) - GRID_MARGIN
    hi = coords.max(0) + GRID_MARGIN
    n = GRID_POINTS
    step = (hi - lo) / (n - 1)
    return lo, step, (n, n, n)


def make_global_grid(all_frames_atoms):
    """Create a universal grid encompassing all frames to prevent resolution changes."""
    all_coords = []
    for atoms in all_frames_atoms:
        all_coords.extend([[x, y, z] for _, x, y, z in atoms])
    coords = np.array(all_coords) * ANG2BOHR
    lo = coords.min(0) - GRID_MARGIN
    hi = coords.max(0) + GRID_MARGIN
    n = GRID_POINTS
    step = (hi - lo) / (n - 1)
    return lo, step, (n, n, n)


def eval_mo_on_grid(mol, mo_vec, origin, step, npts):
    """Evaluate a single MO vector on a 3D grid. Returns flat array."""
    nx, ny, nz = npts
    xs = origin[0] + np.arange(nx) * step[0]
    ys = origin[1] + np.arange(ny) * step[1]
    zs = origin[2] + np.arange(nz) * step[2]
    gx, gy, gz = np.meshgrid(xs, ys, zs, indexing='ij')
    coords = np.column_stack([gx.ravel(), gy.ravel(), gz.ravel()])
    ao = mol.eval_gto('GTOval', coords)
    return ao @ mo_vec


def _cube_header(comment, atoms, origin, step, npts):
    """Generate the header portion of a cube file."""
    nx, ny, nz = npts
    lines = [f'{comment}\n', 'PySCF cube\n']
    lines.append(f'{len(atoms):5d}{origin[0]:12.6f}{origin[1]:12.6f}{origin[2]:12.6f}\n')
    lines.append(f'{nx:5d}{step[0]:12.6f}{0.0:12.6f}{0.0:12.6f}\n')
    lines.append(f'{ny:5d}{0.0:12.6f}{step[1]:12.6f}{0.0:12.6f}\n')
    lines.append(f'{nz:5d}{0.0:12.6f}{0.0:12.6f}{step[2]:12.6f}\n')
    for sym, x, y, z in atoms:
        an = ATOMIC_NUMBERS.get(sym, 0)
        lines.append(f'{an:5d}{float(an):12.6f}{x*ANG2BOHR:12.6f}{y*ANG2BOHR:12.6f}{z*ANG2BOHR:12.6f}\n')
    return lines


def _cube_data(vals, nz):
    """Generate the data portion of a cube file."""
    lines = []
    for i, v in enumerate(vals):
        lines.append(f'{v:13.5e}')
        if i % nz % 6 == 5 or i % nz == nz - 1:
            lines.append('\n')
    return lines


def write_cube_file(fname, comment, atoms, origin, step, npts, vals):
    """Write a Gaussian cube file to disk."""
    with open(fname, 'w') as f:
        f.writelines(_cube_header(comment, atoms, origin, step, npts))
        nz = npts[2]
        idx = 0
        nx, ny, _ = npts
        for ix in range(nx):
            for iy in range(ny):
                for iz in range(nz):
                    f.write(f'{vals[idx]:13.5e}')
                    idx += 1
                    if iz % 6 == 5 or iz == nz - 1:
                        f.write('\n')


def cube_to_string(comment, atoms, origin, step, npts, vals):
    """Generate cube file content as a string (for in-memory compression)."""
    parts = _cube_header(comment, atoms, origin, step, npts)
    nz = npts[2]
    idx = 0
    nx, ny, _ = npts
    for ix in range(nx):
        for iy in range(ny):
            for iz in range(nz):
                parts.append(f'{vals[idx]:13.5e}')
                idx += 1
                if iz % 6 == 5 or iz == nz - 1:
                    parts.append('\n')
    return ''.join(parts)
