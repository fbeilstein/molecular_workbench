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
                isoval = float(max(0.02, np.max(absvals) * 0.12))
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
