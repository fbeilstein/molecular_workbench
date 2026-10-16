# Transition State Engines & Search Strategies

The Transition State Workbench (`server_irc.py`) supports three robust, transparent TS search strategies:

---

## 1. 🧬 Genetic Swarm (Default Intermolecular)

A multi-generation genetic algorithm exploring $SE(3) \times SO(3)$ (both 3D translation and 3D rigid-body quaternions for both reactants).

- **Stage 1 (Global Swarm Exploration)**: Optimizes relative orientation and approach vectors of reactants within an active energy and contact distance window.
- **Stage 2 (Quantum Valley Relaxation)**: Relaxes the top candidates along steepest descent paths toward the reaction valley.
- **Stage 3 (Sella DLC Polish)**: Employs Sella with Delocalized Internal Coordinates (DLC) to locate the exact first-order saddle point (transition state) with a single imaginary frequency.

---

## 2. 🧱 Dual-Endpoint IDPP (Universal All Reaction Classes)

Image-Dependent Pair Potential (IDPP) geometric interpolation between Reactant ($R$) and Product ($P$) endpoints.

- **Universal Interpolation**: Operates on any number of forming and breaking bonds simultaneously (substitutions, Diels-Alder cycloadditions, additions, eliminations, and rearrangements) without reaction-specific heuristics.
- **Atom-Mapping Guided**: Utilizes validated Ketcher 1-to-1 atom indices to establish the reaction coordinate graph.
- **Saddle Convergence**: Initializes a band of intermediate images and relaxes the highest-energy image using Sella saddle optimization into the true transition state.

---

## 3. ⚡ Machine Learning Flow-Matching (RitS / TS-Diff)

Generative neural-network models for predicting transition state geometries directly from reactant and product graphs.

- **Prerequisites Required**: Flow-matching neural network models require PyTorch, PyTorch Geometric, and external model weight checkpoints (~1.5 GB).
- **Zero Silent Failures**: If selected without installed packages or checkpoint files, the engine **never silently fakes or falls back to ground states**. It immediately raises an informative error pointing to setup instructions.
- For complete setup instructions and model weights downloads, see [docs/ML_MODELS.md](ML_MODELS.md).
