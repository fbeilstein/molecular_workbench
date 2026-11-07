# Architecture Overview

Molecular Workbench is an interactive computational chemistry suite for exploring molecular structures, localized orbitals, reaction mechanisms (IRC), and generating lightweight, embeddable 3D visualizations (`.rxnbundle.zip`).

---

## System Components

The suite consists of three specialized server workbenches and a zero-backend standalone web player:

```
┌─────────────────────────────────────────────────────────────────────────────┐
│                            Molecular Workbench                              │
├──────────────────────┬──────────────────────┬───────────────────────────────┤
│ server_multistep.py  │    server_irc.py     │          server.py            │
│ (Multistep Studio)   │ (TS & IRC Workbench) │ (Single Molecule Workbench)   │
│ :8083                │ :8082                │ :4321                         │
├──────────────────────┴──────────────────────┴───────────────────────────────┤
│                              Shared Tools                                   │
│  - tools/qm/ (PySCF DFT, IBO Localization, Topological Network Extraction)   │
│  - irc_plugins/ts_engine/ (Genetic Swarm, Universal IDPP, Sella DLC, RitS)  │
│  - tools/cube_converter.py (Marching Cubes -> 3D Surface JSON Meshes)       │
├─────────────────────────────────────────────────────────────────────────────┤
│                          viewer.html (Standalone)                           │
│  Client-side zero-backend viewer for .rxnbundle.zip archives via WebGL      │
└─────────────────────────────────────────────────────────────────────────────┘
```

---

## 1. Multistep Reaction Pathway Studio (`server_multistep.py`)

A 5-stage interactive studio for building and rendering complex multistep chemical mechanisms:

- **Stage 1 (Trajectory Selection)**: Ingests reaction step JSONs or IRC trajectories.
- **Stage 2 (Stitching & Alignment)**: Aligns multi-fragment coordinate frames across reaction steps using rotational superposition (Kabsch algorithm) and boundary continuity.
- **Stage 3 (Energy Profile Stitching)**: Harmonizes relative electronic energies, barriers, and reaction enthalpies across sequential steps into a single Potential Energy Surface (PES).
- **Stage 4 (Orbital Schedule & ESP)**: Tracks and renders frontier MOs, localized bonds, and continuous Electrostatic Potential (ESP) surfaces across the pathway.
- **Stage 5 (SVG Mechanism Wiring & Bundle Exporter)**: Provides interactive 2D SVG mechanism wiring (reaction arrows, state jumping, electron-pushing arrows) and exports maximal-compression bundles (`.rxnbundle.zip`).

---

## 2. Transition State & IRC Workbench (`server_irc.py`)

A 3-stage workbench for finding reaction transition states and validating reaction paths:

- **Stage 1 (Reaction Input)**: SMILES input or 2D molecular drawing for Reactants and Products.
- **Stage 2 (Atom Mapping)**: Interactive Ketcher 1-to-1 atom index mapping with automated graph validation.
- **Stage 3 (TS Search & IRC)**:
  - Supports Genetic Swarm, Dual-Endpoint IDPP, ORCA, and ML Flow-Matching (RitS).
  - Normal vibrational mode analysis (displacements, single imaginary frequency confirmation).
  - Forward and reverse Intrinsic Reaction Coordinate (IRC) integration.

---

## 3. Single-Molecule Quantum Workbench (`server.py`)

- 3D molecular structure generation and xTB / PySCF geometry optimization.
- Intrinsic Bond Orbital (IBO) localization separating $\sigma$ and $\pi$ systems while preserving textbook Lewis structures.
- Exact virtual subspace antisymmetrization projection for textbook $\sigma^*$ and $\pi^*$ antibonding orbitals.
- Knizia IAO virtual $p_\perp$ projector for clean planar conjugated $\pi^*$ systems free of Rydberg contamination.
- Diagonal Fock operator expectation energies ($E_k = \langle \phi_k | \hat{F} | \phi_k \rangle$, eV) evaluated on all localized bonds, lone pairs, and antibonds.
- Symmetric Löwdin orthonormalization for degenerate canonical eigenspaces ($|\Delta \varepsilon| < 10^{-4}$ Ha) aligned with molecular axes.
- Post-localization overlap-gated snap ($w \ge 0.90$) against canonical frontier MOs (`HOMO-k ... LUMO+k`).
- Unified 0.045 mesh isovalue extraction with multi-isovalue composite bundling (e.g. [18]annulene).
- 3-decimal precision coordinate/mesh trimming, compact JSON packaging, and level-9 compressed `.rxnbundle.zip` export.

---

## 4. Standalone Web Viewer (`viewer.html`)

- **Zero-Backend Architecture**: Runs completely in the browser; can be hosted on GitHub Pages or static web servers.
- **In-Memory Decompression**: Uses `JSZip` to unpack trajectory coordinates, SVG wiring, and precomputed 3D orbital surface meshes directly into memory.
- **Hardware-Accelerated Rendering**: Uses `3Dmol.js` WebGL canvas with synchronized 2D PES canvas playback and interactive mechanism SVG integration.
- **Universal Bundle Support**: Natively loads both multistep reaction bundles (`bundles/reactions/`) and single-molecule orbital archives (`bundles/molecules/`).
