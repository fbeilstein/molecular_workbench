import json
import os
import numpy as np
from skimage import measure
import trimesh

ANG2BOHR = 1.8897259886

class CubeConverter:
    """Processes .cube files AFTER calculation to extract lightweight surface meshes."""
    
    @staticmethod
    def read_cube(filepath):
        with open(filepath, 'r') as f:
            lines = f.readlines()
            
        natoms_line = lines[2].split()
        natoms = int(natoms_line[0])
        origin = np.array([float(x) for x in natoms_line[1:4]])
        
        nx_line = lines[3].split()
        ny_line = lines[4].split()
        nz_line = lines[5].split()
        
        npts = (int(nx_line[0]), int(ny_line[0]), int(nz_line[0]))
        step = np.array([float(nx_line[1]), float(ny_line[2]), float(nz_line[3])])
        
        data_lines = lines[6 + natoms:]
        vals = []
        for line in data_lines:
            vals.extend([float(x) for x in line.split()])
            
        vol = np.array(vals).reshape(npts)
        return vol, origin, step

    @staticmethod
    def extract_mesh(vol, origin, step, level):
        try:
            verts, faces, normals, _ = measure.marching_cubes(vol, level)
            
            # marching_cubes normals point toward higher values.
            # For negative isosurfaces, higher values are outside (e.g. 0 > -0.03), so normals point inwards. We must flip them.
            if level < 0:
                normals = -normals
                faces = faces[:, ::-1]
                
            # Use trimesh to clean up the mesh (remove degenerate faces, unify normals)
            mesh = trimesh.Trimesh(vertices=verts, faces=faces, vertex_normals=normals, process=True)
            
            # Convert back to arrays
            verts = mesh.vertices
            faces = mesh.faces
            normals = mesh.vertex_normals
            
            # Convert to Angstrom coordinates
            verts_angstrom = verts * (step / ANG2BOHR) + (origin / ANG2BOHR)
            
            # Format for 3Dmol.js CustomShape
            v_list = [{'x': float(v[0]), 'y': float(v[1]), 'z': float(v[2])} for v in verts_angstrom]
            f_list = [int(idx) for face in faces for idx in face]
            n_list = [{'x': -float(n[0]), 'y': -float(n[1]), 'z': -float(n[2])} for n in normals]
            
            return {
                'vertices': v_list,
                'faces': f_list,
                'normals': n_list
            }
        except Exception:
            return None

    @classmethod
    def process_cube(cls, cube_path, isoval=0.025, adaptive=False):
        """Converts a .cube file into _pos.json and _neg.json surface files.
        
        If adaptive=True, ignores isoval and instead computes the threshold
        from the 85th percentile of |values|. This keeps orbital volume
        consistent across trajectory frames where normalization fluctuates.
        """
        vol, origin, step = cls.read_cube(cube_path)
        
        if adaptive:
            absvals = np.abs(vol.ravel())
            absvals = absvals[absvals > 1e-8]  # ignore near-zero
            if len(absvals) > 0:
                # Increase minimum isovalue from 0.02 to 0.045 to hide IBO tails visually
                isoval = float(max(0.045, np.max(absvals) * 0.15))
            else:
                isoval = 0.05
        
        extracted = 0
        
        pos_mesh = cls.extract_mesh(vol, origin, step, isoval)
        if pos_mesh:
            pos_path = cube_path.replace('.cube', '_pos.json')
            with open(pos_path, 'w') as f:
                json.dump(pos_mesh, f, separators=(',', ':'))
            extracted += 1
            
        neg_mesh = cls.extract_mesh(vol, origin, step, -isoval)
        if neg_mesh:
            neg_path = cube_path.replace('.cube', '_neg.json')
            with open(neg_path, 'w') as f:
                json.dump(neg_mesh, f, separators=(',', ':'))
            extracted += 1
            
        return extracted > 0

if __name__ == '__main__':
    import sys
    import argparse
    import glob

    parser = argparse.ArgumentParser(description="Interactively select and convert .cube files to .json meshes.")
    parser.add_argument("path", help="Directory containing .cube files, or a specific .cube file")
    parser.add_argument("--adaptive", action="store_true", help="Use adaptive isovalue")
    parser.add_argument("--isoval", type=float, default=0.025, help="Isovalue for mesh extraction")
    args = parser.parse_args()

    if os.path.isdir(args.path):
        cube_files = sorted(glob.glob(os.path.join(args.path, "*.cube")))
    elif args.path.endswith('.cube') and os.path.exists(args.path):
        cube_files = [args.path]
    else:
        print(f"Error: {args.path} is not a valid directory or .cube file.")
        sys.exit(1)

    if not cube_files:
        print("No .cube files found.")
        sys.exit(0)

    print("Found the following orbitals (cube files):")
    for i, cf in enumerate(cube_files, 1):
        print(f"  [{i}] {os.path.basename(cf)}")
    
    print("\nWhich ones do you want to save? (Enter comma-separated numbers, e.g. '1, 3, 4', or 'all')")
    try:
        selection = input("> ").strip()
    except (KeyboardInterrupt, EOFError):
        print("\nCancelled.")
        sys.exit(0)
    
    if selection.lower() == 'all':
        selected_files = cube_files
    else:
        selected_files = []
        try:
            indices = [int(x.strip()) for x in selection.split(',') if x.strip()]
            for idx in indices:
                if 1 <= idx <= len(cube_files):
                    selected_files.append(cube_files[idx - 1])
                else:
                    print(f"Warning: Index {idx} is out of range.")
        except ValueError:
            print("Invalid input. Please enter numbers separated by commas.")
            sys.exit(1)
            
    if not selected_files:
        print("No files selected to process.")
        sys.exit(0)
        
    print(f"\nProcessing {len(selected_files)} files...")
    for cf in selected_files:
        print(f"  -> Converting {os.path.basename(cf)}...")
        CubeConverter.process_cube(cf, isoval=args.isoval, adaptive=args.adaptive)
    print("Done!")
