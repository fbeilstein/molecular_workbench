import numpy as np
from scipy.spatial.transform import Rotation, Slerp
from ase import Atoms
from ase.optimize import LBFGS
from sella import Sella
from .common import rcov, _unit
from .calculators import make_xtb, _make_restrained_calc


def _random_cone_vector(v_axis, max_angle_rad):
    """Sample a random unit vector within a spherical cone around v_axis."""
    z = np.array([0.0, 0.0, 1.0])
    cos_max = np.cos(max_angle_rad)
    cos_theta = np.random.uniform(cos_max, 1.0)
    sin_theta = np.sqrt(max(0.0, 1.0 - cos_theta * cos_theta))
    phi = np.random.uniform(0.0, 2.0 * np.pi)
    v_local = np.array([sin_theta * np.cos(phi), sin_theta * np.sin(phi), cos_theta])
    if np.allclose(v_axis, z):
        return v_local
    rot_align, _ = Rotation.align_vectors([v_axis], [z])
    return _unit(rot_align.apply(v_local))


def _run_unimolecular_rotamer_swarm(atoms, model, log, n_confs=16):
    """Explore rotamer conformations in the reaction valley for unimolecular steps."""
    log("  Unimolecular reaction detected: exploring internal rotamer ensemble in reaction valley...")
    from rdkit import Chem
    from rdkit.Chem import AllChem
    from .calculators import make_xtb, _make_restrained_calc, guard_saddle_polish
    from .vibrations import analyse_mode
    from ase.optimize import LBFGS

    mol_ref = model.r_mols[0] if model.r_mols else None
    confs_pos = []
    if mol_ref is not None:
        try:
            m = Chem.Mol(mol_ref)
            m = Chem.AddHs(m)
            cids = list(AllChem.EmbedMultipleConfs(m, n_confs, AllChem.ETKDGv3()))
            for cid in cids:
                confs_pos.append(m.GetConformer(cid).GetPositions())
        except Exception:
            pass

    if not confs_pos:
        confs_pos = [atoms.positions.copy()]

    targets = [(i, j, model.d0(i, j) + 0.35) for i, j in model.broken]
    for i, j in model.formed:
        targets.append((i, j, model.d0(i, j) + 0.65))

    solvent = getattr(model, "solvent", "water")
    candidates = []

    for idx, p in enumerate(confs_pos):
        at = Atoms(symbols=model.symbols, positions=p)
        base = at.copy()
        base.calc = make_xtb(model.charge, model.multiplicity, solvent=solvent)
        spec = [(i, j, float(np.linalg.norm(p[i] - p[j]))) for i, j in model.spectator_bonds]
        at.calc = _make_restrained_calc(base, targets, k=35.0, spectator_restraints=spec, model=model)
        opt = LBFGS(at, logfile=None, maxstep=0.15)
        try:
            opt.run(fmax=0.08, steps=30)
            base.set_positions(at.get_positions())
            e = float(base.get_potential_energy())
            candidates.append((e, at.get_positions().copy(), idx))
        except Exception:
            pass

    if candidates:
        candidates.sort(key=lambda x: x[0])
        best_e, best_pos, best_idx = candidates[0]
        atoms.set_positions(best_pos)
        de_kcal = (candidates[-1][0] - best_e) * 23.06
        log(f"  ✓ Rotamer swarm evaluated {len(candidates)} conformers: selected lowest-energy valley pose ({best_e:.3f} eV, {de_kcal:.1f} kcal/mol lower).")
    else:
        log("  Proceeding with current geometry.")

    log("  Polishing saddle point with Sella DLC...")
    atoms.calc = make_xtb(model.charge, model.multiplicity, solvent=solvent)
    ts_rep = None
    try:
        ts_rep = analyse_mode(atoms, model, log)
    except Exception:
        pass

    if ts_rep and ts_rep.get("n_imag") == 1 and ts_rep.get("mode_consistent", True):
        log(f"  ✓ Saddle point confirmed by vibrational analysis ({ts_rep['imag_freqs'][0]:.0f}i cm-1).")
        return [{"step": 1, "energy_ev": float(atoms.get_potential_energy()), "fmax": float(np.max(np.linalg.norm(atoms.get_forces(), axis=1)))}]

    return guard_saddle_polish(atoms, model, log, max_steps=20, max_drift=0.15)


def run_genetic_swarm_ts(atoms, model, log, n_gen=8, pop_size=60):
    """Multi-generation Genetic Algorithm exploring SE(3) x SO(3) poses -> Valley Relax -> Sella."""
    log("=== Strategy: Genetic Swarm SE(3) Dual-Rigid-Body Evolution ===")
    if not model.formed or len(model.r_mols) < 2:
        return _run_unimolecular_rotamer_swarm(atoms, model, log)
    solvent = getattr(model, "solvent", "water")

    i_sub, j_nuc = model.formed[0]
    fi_sub = model.frag_of[i_sub]
    fj_nuc = model.frag_of[j_nuc]
    sub_idx = [k for k in range(model.n_atoms) if model.frag_of[k] == fi_sub]
    nuc_idx = [k for k in range(model.n_atoms) if model.frag_of[k] == fj_nuc]
    other_idx = [k for k in range(model.n_atoms) if model.frag_of[k] not in (fi_sub, fj_nuc)]

    p_sub_orig = atoms.positions[i_sub].copy()
    d_ts = 1.45 * (rcov(atoms.symbols[i_sub]) + rcov(atoms.symbols[j_nuc]))

    v_user = _unit(atoms.positions[j_nuc] - p_sub_orig)
    if np.linalg.norm(v_user) < 1e-4:
        v_user = np.array([1.0, 0.0, 0.0])

    p_sub_local = atoms.positions[sub_idx] - atoms.positions[i_sub]
    p_nuc_local = atoms.positions[nuc_idx] - atoms.positions[j_nuc]

    def decode(gene):
        r_val, u_vec, q_sub, q_nuc = gene
        r_mat_sub = Rotation.from_quat(q_sub)
        r_mat_nuc = Rotation.from_quat(q_nuc)

        pos = np.zeros((model.n_atoms, 3))
        pos[sub_idx] = r_mat_sub.apply(p_sub_local) + p_sub_orig
        pos[nuc_idx] = r_mat_nuc.apply(p_nuc_local) + p_sub_orig + r_val * r_mat_sub.apply(u_vec)
        for k in other_idx:
            pos[k] = atoms.positions[k]
        return pos

    def evaluate(pos):
        for a in range(model.n_atoms):
            for b in range(a + 1, model.n_atoms):
                if (a, b) in model.formed or (b, a) in model.formed:
                    continue
                if (a, b) in model.r_bonds or (b, a) in model.r_bonds:
                    continue
                d = np.linalg.norm(pos[a] - pos[b])
                if d < 0.75 * (rcov(atoms.symbols[a]) + rcov(atoms.symbols[b])):
                    return 1e6
        test = atoms.copy()
        test.set_positions(pos)
        test.calc = make_xtb(model.charge, model.multiplicity, solvent=solvent)
        try:
            e_base = float(test.get_potential_energy())
        except Exception:
            return 1e6

        e_pen = 0.0
        # Multi-bond guidance: penalize deviation for all forming bonds
        for a, b in model.formed:
            if (a in sub_idx and b in nuc_idx) or (b in sub_idx and a in nuc_idx):
                d_ab = np.linalg.norm(pos[a] - pos[b])
                d_tgt = 1.45 * (rcov(atoms.symbols[a]) + rcov(atoms.symbols[b]))
                e_pen += 15.0 * (d_ab - d_tgt) ** 2
        return e_base + e_pen

    pop = []
    for idx in range(pop_size):
        r_init = float(np.random.uniform(d_ts - 0.20, d_ts + 0.35))
        if idx < int(0.50 * pop_size):
            u_init = _random_cone_vector(v_user, np.radians(60))
        else:
            u_init = _unit(np.random.normal(size=3))

        # Substrate rotation: 50% random, 50% small tilt
        if idx % 2 == 0:
            q_sub = Rotation.random().as_quat()
        else:
            q_sub = Rotation.from_rotvec(np.random.normal(0, np.radians(25), size=3)).as_quat()

        # Nucleophile rotation: 50% uniform Haar, 25% azimuthal sweep, 25% wide tilt
        if idx % 4 in (0, 1):
            q_nuc = Rotation.random().as_quat()
        elif idx % 4 == 2:
            azimuth = (idx * np.pi / 3.0) % (2.0 * np.pi)
            q_nuc = Rotation.from_rotvec(azimuth * u_init).as_quat()
        else:
            q_nuc = Rotation.from_rotvec(np.random.normal(0, np.radians(60), size=3)).as_quat()

        pop.append([r_init, u_init, q_sub, q_nuc])

    log(f"  SE(3) Population: {pop_size} individuals exploring rotations & translations in {solvent.upper()}.")
    for gen in range(n_gen):
        scores = [evaluate(decode(ind)) for ind in pop]
        ranked = np.argsort(scores)
        best_e = scores[ranked[0]]
        best_gene = pop[ranked[0]]
        rot_a = np.degrees(np.linalg.norm(Rotation.from_quat(best_gene[2]).as_rotvec()))
        rot_b = np.degrees(np.linalg.norm(Rotation.from_quat(best_gene[3]).as_rotvec()))
        log(f"  Cycle {gen + 1}/{n_gen}: Best E = {best_e:.3f} eV, r = {best_gene[0]:.2f} A, rotA = {rot_a:.1f}°, rotB = {rot_b:.1f}°")

        new_pop = [pop[k] for k in ranked[:6]]
        while len(new_pop) < pop_size:
            p1 = pop[np.random.choice(ranked[:15])]
            p2 = pop[np.random.choice(ranked[:15])]
            r_c = float(0.5 * (p1[0] + p2[0]) + np.random.normal(0, 0.04))
            u_c = _unit(0.5 * (p1[1] + p2[1]) + np.random.normal(0, 0.08, size=3))
            try:
                slerp_sub = Slerp([0.0, 1.0], Rotation.from_quat([p1[2], p2[2]]))
                q_sub_c = slerp_sub(np.random.uniform(0.2, 0.8)).as_quat()
            except Exception:
                q_sub_c = p1[2].copy()
            try:
                slerp_nuc = Slerp([0.0, 1.0], Rotation.from_quat([p1[3], p2[3]]))
                q_nuc_c = slerp_nuc(np.random.uniform(0.2, 0.8)).as_quat()
            except Exception:
                q_nuc_c = p1[3].copy()

            # Multi-scale mutation for both fragments
            w_sub = np.random.normal(0, np.radians(45 if np.random.rand() < 0.25 else 12), size=3)
            q_sub_m = (Rotation.from_rotvec(w_sub) * Rotation.from_quat(q_sub_c)).as_quat()
            w_nuc = np.random.normal(0, np.radians(60 if np.random.rand() < 0.35 else 15), size=3)
            q_nuc_m = (Rotation.from_rotvec(w_nuc) * Rotation.from_quat(q_nuc_c)).as_quat()
            new_pop.append([r_c, u_c, q_sub_m, q_nuc_m])
        pop = new_pop

    scores = [evaluate(decode(ind)) for ind in pop]
    best_idx = int(np.argmin(scores))
    atoms.set_positions(decode(pop[best_idx]))
    d_final = float(np.linalg.norm(atoms.positions[i_sub] - atoms.positions[j_nuc]))
    log(f"  Swarm converged to lowest-strain 3D pose ({d_final:.3f} A).")

    # Stage 2: Quantum Valley Internal Relaxation (relaxing umbrella & dihedrals)
    log("  Stage 2: Quantum valley relaxation of internal coordinates...")
    base = atoms.copy()
    base.calc = make_xtb(model.charge, model.multiplicity, solvent=solvent)
    targets = [(i, j, model.d0(i, j) + 0.65) for i, j in model.formed]
    targets += [(i, j, model.d0(i, j) + 0.35) for i, j in model.broken]
    pos_cur = atoms.positions.copy()
    spec = [(i, j, float(np.linalg.norm(pos_cur[i] - pos_cur[j]))) for i, j in model.spectator_bonds]
    calc = _make_restrained_calc(base, targets, k=35.0, spectator_restraints=spec, model=model)
    atoms.calc = calc
    opt_valley = LBFGS(atoms, logfile=None, maxstep=0.15)
    opt_valley.run(fmax=0.08, steps=40)
    log(f"  Internal coordinates relaxed in {opt_valley.nsteps} steps into quadratic valley.")

    # Stage 3: Sella DLC Saddle Point Polish (guarded)
    log("  Stage 3: Polishing saddle point with Sella in internal coordinates...")
    from .calculators import guard_saddle_polish
    from .vibrations import analyse_mode
    atoms.calc = make_xtb(model.charge, model.multiplicity, solvent=solvent)
    ts_rep = None
    try:
        ts_rep = analyse_mode(atoms, model, log)
    except Exception:
        pass

    if ts_rep and ts_rep.get("n_imag") == 1 and ts_rep.get("mode_consistent", True):
        log(f"  ✓ Saddle point confirmed by vibrational analysis ({ts_rep['imag_freqs'][0]:.0f}i cm-1).")
        return [{"step": 1, "energy_ev": float(atoms.get_potential_energy()), "fmax": float(np.max(np.linalg.norm(atoms.get_forces(), axis=1)))}]

    return guard_saddle_polish(atoms, model, log, max_steps=50, max_drift=0.35)

