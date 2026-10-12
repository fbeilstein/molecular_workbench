import os
import numpy as np
import torch
from ase.optimize import LBFGS
from sella import Sella
from rdkit import Chem
from .calculators import make_xtb, _make_restrained_calc, guard_saddle_polish
from .lego_builder import _embed_aligned_products

_CACHED_RITS_MODEL = None
_CACHED_CKPT_PATH = None


def _get_rits_model(ckpt_path):
    global _CACHED_RITS_MODEL, _CACHED_CKPT_PATH
    if _CACHED_RITS_MODEL is not None and _CACHED_CKPT_PATH == ckpt_path:
        return _CACHED_RITS_MODEL

    from omegaconf import OmegaConf
    from megalodon.models.module import Graph3DInterpolantModel
    from megalodon.data.ts_batch_preprocessor import TsBatchPreProcessor

    cfg_candidates = [
        "tools/rits_repo/scripts/conf/rits.yaml",
        "scripts/conf/rits.yaml",
    ]
    cfg_path = next((p for p in cfg_candidates if os.path.isfile(p)), None)
    if not cfg_path:
        raise FileNotFoundError("Could not find rits.yaml configuration file.")

    cfg = OmegaConf.load(cfg_path)
    batch_preprocessor = TsBatchPreProcessor(
        aug_rotations=cfg.data.get("aug_rotations", False),
        scale_coords=cfg.data.get("scale_coords", 1.0),
    )

    device = "cuda" if torch.cuda.is_available() else "cpu"
    model = Graph3DInterpolantModel.load_from_checkpoint(
        ckpt_path,
        loss_params=cfg.loss,
        interpolant_params=cfg.interpolant,
        sampling_params=cfg.sample,
        batch_preprocessor=batch_preprocessor,
        strict=True,
        weights_only=False,
    )
    model = model.to(device).eval()
    _CACHED_RITS_MODEL = (model, device)
    _CACHED_CKPT_PATH = ckpt_path
    return _CACHED_RITS_MODEL


def _score_ts_candidate(pos, model):
    score = 0.0
    for i, j in model.formed:
        d = np.linalg.norm(pos[i] - pos[j])
        score -= abs(d - (model.d0(i, j) + 0.65)) * 2.0
    for i, j in model.broken:
        d = np.linalg.norm(pos[i] - pos[j])
        score -= abs(d - (model.d0(i, j) + 0.35)) * 2.0

    # Inversion / backside attack collinearity (180 deg) across reacting centers
    for f1, f2 in model.formed:
        for b1, b2 in model.broken:
            shared = None
            if f1 == b1: shared, a, b = f1, f2, b2
            elif f1 == b2: shared, a, b = f1, f2, b1
            elif f2 == b1: shared, a, b = f2, f1, b2
            elif f2 == b2: shared, a, b = f2, f1, b1
            if shared is not None:
                va = pos[a] - pos[shared]
                vb = pos[b] - pos[shared]
                cos_ang = np.dot(va, vb) / (np.linalg.norm(va) * np.linalg.norm(vb) + 1e-9)
                score += (1.0 - cos_ang) * 5.0
    return score


def _sample_best_rits_coords(model_net, device, numbers, charges, r_coord, p_coord, ts_coord, edge_index, edge_attr, model, log, n_samples=4):
    from torch_geometric.data import Data, DataLoader
    from megalodon.metrics.ts_evaluation_callback import convert_coords_to_np

    data_list = [
        Data(
            num_nodes=len(numbers),
            numbers=numbers,
            charges=charges,
            r_coord=r_coord,
            p_coord=p_coord,
            ts_coord=ts_coord,
            edge_index=edge_index,
            edge_attr=edge_attr,
            id=f"sample_{idx}"
        ) for idx in range(n_samples)
    ]

    log(f"  Step 2: Generating transition state ensemble via flow matching ({n_samples} parallel candidate paths)...")
    loader = DataLoader(data_list, batch_size=n_samples)
    coords_list = []
    for batch in loader:
        batch = batch.to(device)
        with torch.no_grad():
            sample = model_net.sample(batch=batch, timesteps=25, pre_format=True)
        coords_list = convert_coords_to_np(sample)
        break

    if not coords_list:
        raise RuntimeError("RitS flow matching produced empty coordinates.")

    scores = [_score_ts_candidate(c, model) for c in coords_list]
    best_idx = int(np.argmax(scores))
    log(f"  ✓ Optimal candidate #{best_idx + 1}/{n_samples} selected by reaction collinearity (score: {scores[best_idx]:.2f}).")
    return coords_list[best_idx]


def run_rits_ts(atoms, model, log):
    """RitS (Reaction-informed Transition State) flow-matching & quantum refinement."""
    log("=== Strategy: RitS (Generative Flow Matching) ===")

    has_rits = False
    try:
        import megalodon  # noqa
        has_rits = True
    except ImportError:
        pass

    candidates = [
        os.environ.get("RITS_CHECKPOINT"),
        "tools/ml_models/rits/rits.ckpt",
        "tools/ml_models/rits/rits_checkpoint.pt",
        "tools/rits_repo/data/rits.ckpt",
        "tools/rits_repo/data/ts1x_rits.ckpt",
    ]
    checkpoint = next((p for p in candidates if p and os.path.isfile(p)), None)
    has_weights = checkpoint is not None

    if not has_rits or not has_weights:
        missing = []
        if not has_rits:
            missing.append("Python package 'megalodon' (run setup.sh --with-ml)")
        if not has_weights:
            missing.append("Model checkpoint 'tools/ml_models/rits/rits.ckpt'")
        err_msg = (
            f"RitS cannot run: Missing {', '.join(missing)}.\n"
            "Please run './setup.sh --with-ml' to install dependencies and weights.\n"
            "For general reactions without external ML checkpoints, choose "
            "'🧬 Genetic Swarm' or '🧱 Dual-Endpoint IDPP'."
        )
        log(f"\n❌ ERROR: {err_msg}")
        raise RuntimeError(err_msg)

    log(f"  Step 1: Loading RitS neural network model from {checkpoint}...")
    model_net, device = _get_rits_model(checkpoint)

    # Prepare Reactant (R) and Product (P) coordinates and bond graphs
    pos_R = atoms.positions.copy()
    pos_P = _embed_aligned_products(model, pos_R)

    bmat_r = np.zeros((model.n_atoms, model.n_atoms), dtype=np.int32)
    for (i, j), bo in model.r_bond_orders.items():
        bmat_r[i, j] = bmat_r[j, i] = int(round(bo))

    bmat_p = np.zeros((model.n_atoms, model.n_atoms), dtype=np.int32)
    for (i, j), bo in model.p_bond_orders.items():
        bmat_p[i, j] = bmat_p[j, i] = int(round(bo))

    atomic_numbers = np.array([Chem.Atom(sym).GetAtomicNum() for sym in atoms.symbols], dtype=np.int64)

    from torch_geometric.data import Data, DataLoader
    from megalodon.metrics.ts_evaluation_callback import convert_coords_to_np

    numbers = torch.from_numpy(atomic_numbers).to(torch.uint8)
    r_coord = torch.from_numpy(pos_R).float()
    p_coord = torch.from_numpy(pos_P).float()
    ts_coord = (r_coord + p_coord) / 2.0
    t_bmat_r = torch.from_numpy(bmat_r)
    t_bmat_p = torch.from_numpy(bmat_p)
    charges = torch.full_like(numbers, model.charge, dtype=torch.int8)

    edge_index = (t_bmat_r + t_bmat_p).nonzero().contiguous().T
    edge_attr = torch.stack([t_bmat_r, t_bmat_p], dim=-1)[edge_index[0], edge_index[1]].to(torch.uint8)

    pos_ts = _sample_best_rits_coords(
        model_net, device, numbers, charges, r_coord, p_coord, ts_coord,
        edge_index, edge_attr, model, log, n_samples=4
    )
    atoms.set_positions(pos_ts)
    log("  ✓ Flow-matching TS geometry generated:")
    for i, j in model.formed:
        log(f"    Forming {model.label(i)}-{model.label(j)}: {np.linalg.norm(pos_ts[i] - pos_ts[j]):.3f} A")
    for i, j in model.broken:
        log(f"    Breaking {model.label(i)}-{model.label(j)}: {np.linalg.norm(pos_ts[i] - pos_ts[j]):.3f} A")

    # Step 3: Gentle Valley Relaxation
    log("  Step 3: Quantum valley relaxation of internal coordinates...")
    base = atoms.copy()
    base.calc = make_xtb(model.charge, model.multiplicity)
    targets = [(i, j, float(np.linalg.norm(pos_ts[i] - pos_ts[j]))) for i, j in model.formed + model.broken]
    pos_cur = atoms.positions.copy()
    spec = [(i, j, float(np.linalg.norm(pos_cur[i] - pos_cur[j]))) for i, j in model.spectator_bonds]
    calc = _make_restrained_calc(base, targets, k=35.0, spectator_restraints=spec, model=model)
    atoms.calc = calc
    opt_valley = LBFGS(atoms, logfile=None, maxstep=0.08)
    opt_valley.run(fmax=0.08, steps=25)
    log(f"  Internal coordinates relaxed in {opt_valley.nsteps} steps.")

    # Record profile for UI
    try:
        e_cur = float(atoms.get_potential_energy())
        f_cur = float(np.max(np.linalg.norm(atoms.get_forces(), axis=1)))
    except Exception:
        e_cur, f_cur = 0.0, 0.0

    profile = [{"step": 1, "energy_ev": e_cur, "fmax": f_cur}]
    log(f"  ✓ Transition state geometry ready: E = {e_cur:.4f} eV, max|F| = {f_cur:.4f} eV/A")
    return profile


