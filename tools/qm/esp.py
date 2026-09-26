import os
import numpy as np
import json
from pyscf import gto, tools
from skimage import measure
import trimesh

ANG2BOHR = 1.8897259886


def compute_esp_surface(mol, mf, name, out_dir, grid_points=50, iso_density=0.01):
    """Compute molecular electrostatic potential mapped onto electron density isosurface.

    Steps:
      1. Compute total electron density on a 3D grid
      2. Extract isosurface at iso_density using marching cubes
      3. Evaluate the electrostatic potential at each surface vertex
      4. Save mesh + per-vertex ESP values as compact JSON

    The frontend can then color the surface using a red↔blue diverging map.
    """
    from pyscf import gto, tools
    from pyscf.tools import cubegen
    from skimage import measure
    import trimesh

    print(f"  Computing ESP surface (isodensity={iso_density})...")

    # Build the density matrix
    dm = mf.make_rdm1()
    if isinstance(dm, (tuple, list)):
        dm = dm[0] + dm[1]  # alpha + beta for open-shell

    # Generate electron density cube
    dens_cube = os.path.join(out_dir, f"{name}_density.cube")
    cubegen.density(mol, dens_cube, dm,
                          nx=grid_points, ny=grid_points, nz=grid_points, margin=5.0)

    # Read the cube file
    with open(dens_cube, 'r') as f:
        lines = f.readlines()

    natoms_line = lines[2].split()
    natoms = int(natoms_line[0])
    origin = np.array([float(x) for x in natoms_line[1:4]])

    nx_line = lines[3].split()
    ny_line = lines[4].split()
    nz_line = lines[5].split()
    npts = (int(nx_line[0]), int(ny_line[0]), int(nz_line[0]))
    step = np.array([float(nx_line[1]), float(ny_line[2]), float(nz_line[3])])

    vals = []
    for line in lines[6 + natoms:]:
        vals.extend([float(x) for x in line.split()])
    vol = np.array(vals).reshape(npts)

    # Clean up cube file
    os.remove(dens_cube)

    # Extract isosurface
    try:
        verts, faces, normals, _ = measure.marching_cubes(vol, iso_density)
    except Exception as e:
        print(f"  ✗ ESP surface extraction failed: {e}")
        return None

    # Clean with trimesh
    mesh = trimesh.Trimesh(vertices=verts, faces=faces, vertex_normals=normals, process=True)
    verts = mesh.vertices
    faces = mesh.faces
    normals = mesh.vertex_normals

    # Convert from grid indices to Bohr, then to Angstrom
    verts_bohr = verts * step + origin
    verts_ang = verts_bohr / ANG2BOHR

    # Evaluate ESP at each surface vertex (in Bohr coordinates)
    print(f"  Evaluating ESP at {len(verts)} surface vertices...")

    # Nuclear contribution (vectorized)
    nuc_charges = mol.atom_charges()
    nuc_coords = mol.atom_coords()  # Bohr

    esp_values = np.zeros(len(verts_bohr))
    for a in range(mol.natm):
        dists = np.linalg.norm(verts_bohr - nuc_coords[a], axis=1)
        dists = np.maximum(dists, 1e-6)  # avoid division by zero
        esp_values += nuc_charges[a] / dists

    # Electronic contribution: V_elec(r) = -Tr[D * (1/|r-R|)]
    # Using PySCF's with_rinv_origin + int1e_rinv integral
    print(f"  Evaluating electronic ESP...")
    for i, coord in enumerate(verts_bohr):
        with mol.with_rinv_origin(coord):
            ints = mol.intor('int1e_rinv')
        esp_values[i] -= np.einsum('ij,ij->', dm, ints)

    # Convert ESP from Hartree/e to kcal/mol for chemical intuition
    HARTREE2KCAL = 627.509
    esp_kcal = esp_values * HARTREE2KCAL

    # Build compact JSON mesh
    v_list = [{'x': round(float(v[0]), 4), 'y': round(float(v[1]), 4), 'z': round(float(v[2]), 4)}
              for v in verts_ang]
    f_list = [int(idx) for face in faces for idx in face]
    n_list = [{'x': -round(float(n[0]), 4), 'y': -round(float(n[1]), 4), 'z': -round(float(n[2]), 4)}
              for n in normals]
    esp_list = [round(float(v), 2) for v in esp_kcal]

    esp_data = {
        'vertices': v_list,
        'faces': f_list,
        'normals': n_list,
        'esp_values': esp_list,  # kcal/mol per vertex
        'esp_min': round(float(esp_kcal.min()), 2),
        'esp_max': round(float(esp_kcal.max()), 2),
    }

    esp_path = os.path.join(out_dir, f"{name}_esp.json")
    with open(esp_path, 'w') as f:
        json.dump(esp_data, f, separators=(',', ':'))

    sz = os.path.getsize(esp_path) / 1024
    print(f"  ✓ ESP surface: {len(verts)} vertices, ESP range [{esp_kcal.min():.1f}, {esp_kcal.max():.1f}] kcal/mol → {name}_esp.json ({sz:.0f} KB)")

    return {'file': f"{name}_esp.json", 'esp_min': esp_data['esp_min'], 'esp_max': esp_data['esp_max']}


