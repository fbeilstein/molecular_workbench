# Machine Learning Reaction Path Development Plan

This document outlines the architecture and phased development plan for upgrading the Molecular Workbench to support automated, robust Intrinsic Reaction Coordinate (IRC) mapping using modern Machine Learning techniques.

## Architectural Context
The old approach (manual RDKit atom-mapping + arbitrary driving coordinates + xTB relaxed scans) was brittle and resulted in fragmented/exploding molecules. 
The new architecture relies on state-of-the-art Deep Learning models to handle mapping and 3D positioning, and rigorous optimizers for the physical pathway.

**The Pipeline (one elementary step):**
`Ketcher (R → P, charge, spin)` -> `rxnmapper / Ketcher map + H completion` (2D Map) -> `RitS` (3D TS Guess) **or** fallback `xTB R/P opt → IDPP → CI-NEB` -> `Sella + MLIP/xTB` (TS Optimization) -> `Sella IRC` -> `PySCF IBO tracking` (Orbitals) -> `Electron-flow arrows`

**The Pipeline (multi-step mechanism):** an ordered chain of the above, glued at shared intermediates (see Phase 7).

---

## Feasibility Review (2026-10-04)

**Verdict:** Realistic for **single, closed-shell, organic elementary steps of ≤ ~30–40 atoms**. It is **not** an "input overall reaction → get mechanism" machine. Multi-step chemistry (Grignard, acetal formation, ester hydrolysis) works only if the user supplies the intermediates, with each step computed separately and then stitched together.

| Phase | Original choice | Realism | Main issue / correction |
|---|---|---|---|
| 1. Mapping | `rxnmapper` | ✅ High | Maps heavy atoms only. H mapping must be derived separately. Trained on overall patent reactions, not elementary steps. Ketcher's manual mapping should be the override. |
| 2. TS guess | `OA-ReactDiff` | ⚠️ Low | Trained on Transition1x: **C/H/N/O only, ≤7 heavy atoms, neutral**. Needs **3D** R and P, not just 2D. Replace with **RitS** (graph-only input, CHNOSFP, charged OK) and add a physics fallback (NEB). |
| 3. TS opt | `Sella` + `MACE`/`TorchANI` | ✅ Medium | MACE-OFF / ANI-2x are neutral-organic only (no Mg, Li, B, Al, metals) and not trained on TS regions. Need `GFN2-xTB` as a universal fallback calculator (already in the project). |
| 4. IRC | "ASE IRC module" | ✅ High | **ASE has no IRC.** Use `sella.IRC`. IRC ends are not optimized minima, so add end-point optimization plus a graph check. |
| 5. Orbitals | PySCF, Boys/PM | ✅ Medium (slow) | Should use **IBO** (the project's established choice; see `docs/scientific_choices.md` §4). Cost is minutes to hours per step, so results must be precomputed, not real-time. Classification heuristics will flicker near the TS. |
| — | (missing) | — | No concept of multi-step mechanisms, charge/spin input, solvent, barrierless steps, or electron-flow (curved-arrow) extraction. Added below as Phases 6–7. |

---

## Scope: How Complex Can an Input Reaction Be?

| Tier | What | Route | Examples |
|---|---|---|---|
| **A: automatic** | 1 elementary step, closed-shell, CHNOSFP, neutral or charged, ≤ ~40 atoms | ML route (RitS → MLIP/xTB) | Diels–Alder (endo/exo), Claisen/Cope, SN2, E2, ene, [1,5]-H shift, HCN→HNC |
| **B: semi-automatic** | 1 step with main-group metals (Li, Mg, B, Al), Br/I, a few explicit solvent molecules, ≤ ~80 atoms | Physics fallback (xTB + NEB + Sella) | Grignard C–C formation, organolithium addition, hydroboration, Wittig |
| **C: user-composed** | Multi-step mechanism, every intermediate drawn by the user | Chain of A/B steps + Phase 7 gluing | Full Grignard sequence, acetal formation, ester hydrolysis, aldol condensation, SN1 |
| **D: out of scope** | Open-shell transition metals, radical chains, photochemistry, surface reactions, enzymes, dynamic solvent participation | — | Pd cross-coupling, Grignard *formation* on Mg metal, SET mechanisms, proton shuttles in bulk water |

**Hard limits come from Phase 5 cost:** B3LYP/6-31G* scales steeply with system size. At ~30 atoms a frame takes roughly 0.5–2 min, and an IRC needs 40–100 frames. Anything much above ~40 atoms becomes an overnight job.

---

## Core Concept: The Elementary Step Is the Unit of Work
- An IRC connects **one TS** to **two minima**. The QM machinery never sees "the reaction", only one step at a time.
- A 2D overall equation (e.g., `RMgBr + ketone → alcohol`) does **not** define a mechanism. The user must draw each elementary step (R → P with exactly one TS), or explicitly mark a step as barrierless.
- Barrierless steps (coordination/association, ion recombination, many proton transfers in solution) have **no saddle point**. The pipeline must detect "no TS found" and fall back to a relaxed distance scan for animation instead of failing.

---

## Phase 1: Robust 2D Atom Mapping
- **Goal**: Guarantee flawless atom-to-atom matching between reactants and products, even for complex rearrangements or migrating hydrogens.
- **Tool**: `rxnmapper` (IBM RXN transformer model) or similar ML-driven mapping.
- **Review notes**:
  - `rxnmapper` maps **heavy atoms only**. Migrating hydrogens (tautomerism, H-shifts, proton transfer) require an explicit **H-completion step**: assign H indices from the change in heavy-atom H-counts, and pick the H whose geometry places it closest to its acceptor.
  - It was trained on USPTO *overall* reactions. On user-drawn elementary steps with charges and metals it may be off-distribution. **Ketcher's built-in atom-mapping tool is the manual override** (user-in-the-loop point #2).
  - Charge and spin multiplicity must be entered per step (UI fields) and carried through all phases.
  - For multi-step mechanisms, every step must share **one global atom index space** (see Phase 7).
- **Testing Method**: 
  - **Visual Verification Output**: The script will take the mapped SMILES and generate a 2D image using RDKit (`Draw.ReactionToImage`). The atoms in both the reactant and product will be explicitly labeled with matching integer indices.
  - **Acceptance Criteria**: The user can draw a reaction in Ketcher, hit "Test Mapping", and see an image where atom `1` on the left is undeniably the exact same atom `1` on the right. H atoms involved in bond changes are also labeled.

## Phase 2: Direct 3D Transition State Prediction
- **Goal**: Predict the 3D geometry of the Transition State directly from the mapped 2D structure, bypassing the need to mathematically position reactants in a pre-reaction complex.
- **Tool (primary)**: `RitS` (flow matching, input = reactant/product **connectivity only**, CHNOSFP, neutral + charged, stereo-aware endo/exo and E/Z). GPU recommended (8–24 GB VRAM). CPU is possible but not studied by the authors.
- **Tool (secondary)**: `OA-ReactDiff` / `React-OT`. Only for C/H/N/O, ≤7 heavy atoms, and they need 3D R and P geometries, so they **do not** bypass reactant positioning.
- **Tool (fallback, Phase 2b, any element)**: build 3D R and P complexes with the mapped atom order, optimize both with `GFN2-xTB`, interpolate with **IDPP**, then run **CI-NEB** (ASE `NEB` + xTB). The climbing image becomes the TS guess. This route handles Mg/Li/B/Al and anything xTB covers, at the cost of minutes instead of seconds.
- **Review notes**:
  - Generative models are stochastic. Sample **N = 8–16** TS guesses, optimize all in Phase 3, cluster them, and keep the lowest barrier. Show alternatives to the user if they differ chemically (user-in-the-loop point #3).
  - TS quality depends on reactant conformer. Endo/exo and E/Z must be specified via stereo in the drawn product.
- **Testing Method**: 
  - **3D Render Verification**: The script will output the predicted TS as a standalone `.xyz` file.
  - **Acceptance Criteria**: The TS `.xyz` will be rendered in the 3Dmol.js viewer. The user will visually confirm that the distances and angles of the breaking/forming bonds look like a plausible transition state geometry.

## Phase 3: Rigorous TS Optimization, Verification & Manual Fallback
- **Goal**: Take the AI-predicted TS and optimize it to a true First-Order Saddle Point. If the automated AI route fails or is rejected, provide a 3D editor to manually construct the pre-reaction complex or TS guess.
- **Tools**: 
  - `MACE`, `TorchANI`, or `GFN2-xTB` (for energy/forces).
  - `Sella` (Saddle-point optimizer).
  - **Manual Fallback UI**: A 3D manipulation UI (e.g., `3Dmol.js` or `NGLViewer`) with tools to:
    1. Generate 3D conformers from 2D (e.g., RDKit ETKDG).
    2. Perform quick individual molecule optimizations (e.g., MMFF94 or GFN2-xTB).
    3. Rigid-body rotation/translation of independent fragments to dock them.
    4. **Directly construct the TS Guess**: Manually arrange atoms (e.g., forming a 6-membered ring) and optionally set constrained distances (e.g., "freeze this C-C distance at 1.9 Å") to build the TS by hand, bypassing NEB completely.
    5. Add explicit molecules (solvent, second equivalents).
- **Review notes**:
  - **Calculator coverage**: MACE-OFF (H C N O F P S Cl Br I) and ANI-2x are organic-only. `GFN2-xTB` is the universal fallback for metals (like Mg in Grignard).
  - **The Grignard Edge-Case (Fallback)**: As seen in complex multi-molecular transition states (e.g., 6-membered ring Grignard TS with 2 equivalents of R-MgX), generative AI often fails. If the user rejects the AI output, they enter the 3D editor to manually arrange the 3D fragments into the correct orbital attack pose and re-submit for NEB/optimization.
  - Energies along the path should be reported at **one consistent level**.
  - Frequency check tolerance: small spurious imaginary modes (< ~50i cm⁻¹) are treated as noise.
- **Testing Method**: 
  - **Hessian Frequency Check**: MUST yield exactly ONE significant imaginary frequency.
  - **Acceptance Criteria**: The user watches the vibration in the UI. If the atoms are vibrating along the exact bond breaking/forming coordinate, the TS is approved. If not, the user clicks "Reject" and drops into the 3D manual assembly editor to fix the pose and try again.


## Phase 4: Intrinsic Reaction Coordinate (IRC) Generation
- **Goal**: Slide downhill from the verified TS to generate the full 0K reaction trajectory connecting the Reactant Complex and Product Complex.
- **Tools**: `sella.IRC` (ASE itself has **no** IRC implementation) + the Phase 3 calculator.
- **Review notes**:
  - IRC endpoints are not fully converged minima. **Optimize both ends**, then run an **automatic graph check**: the connectivity of each end must match the drawn R and P (via RDKit `DetermineBonds` / distance criteria). A mismatch means the TS belongs to a different reaction, so reject it and try the next candidate from Phase 2.
  - IRC starts from the *pre-reactive complex*, not from separated molecules. For the "reactants approaching" visual, prepend a short relaxed approach scan from separated fragments.
  - IRC steps are uneven in mass-weighted arc length. **Resample** to uniform spacing for smooth scrubbing.
- **Testing Method**: 
  - **Trajectory Scrubbing**: The script generates a multi-frame `.xyz` trajectory.
  - **Acceptance Criteria**: The user loads the trajectory in the UI and scrubs a timeline slider. They must see a smooth, physical animation of reactants colliding and transforming into products with absolutely zero "exploding" atoms. The automated graph check passes on both ends.

## Phase 5: Continuous Orbital Tracking
- **Goal**: Compute wavefunctions across the entire IRC trajectory and extract specific orbitals (e.g., tracking the HOMO or a specific $\sigma$ bond) as they morph.
- **Tool**: `PySCF` with **Intrinsic Bond Orbitals (IBO)**, consistent with the rest of the pipeline. (The original draft said Boys/Pipek-Mezey. Boys produces banana bonds, which `scientific_choices.md` §4 explicitly rejects.)
- **Review notes**:
  - **Frame-to-frame identity**: match IBOs between consecutive frames by maximum overlap $|\langle\phi_i^{(k)}|S|\phi_j^{(k+1)}\rangle|$ (Hungarian assignment) and fix phase signs. Do not rely on energetic ordering.
  - **Speed**: seed each SCF with the previous frame's density. IBOs are nearly basis-insensitive (they project onto a minimal basis), so a smaller basis (def2-SVP) is acceptable for tracking. Everything is computed **offline** and cached in the bundle.
  - **Classification flicker**: the 10% bond/lone-pair threshold, the σ/π circle test and the hypervalent demotion all assume Lewis-like geometries. Near the TS, orbitals are legitimately "in transit". Assign labels at the **endpoints**, carry identity by overlap, and tag mid-path orbitals as `transforming` instead of re-classifying every frame. (New `.type` strings must be added through `renderer.py` and `bundle_exporter.py`; see `agent_knowledge.md` §4.)
  - σ*/π* projection is ill-defined mid-reaction. Track **occupied** IBOs only along the path.
  - Open-shell steps need UKS and separate α/β IBOs (fishhook arrows). Defer this until after closed-shell works.
- **Testing Method**: 
  - **Volumetric Animation**: Generate surface meshes (`.json`) for the tracked orbital at each frame.
  - **Acceptance Criteria**: As the user scrubs the 3D timeline, the orbital surface dynamically updates without suddenly vanishing or swapping to a completely unrelated part of the molecule.

## Phase 6 (New): Electron Flow — Curved Arrows from First Principles
- **Goal**: Derive textbook curved arrows automatically from the computed wavefunction, with no hand-drawn arrows required.
- **Method**: Knizia & Klein, *"Electron Flow in Reaction Mechanisms—Revealed from First Principles"*, Angew. Chem. Int. Ed. 2015. Along an IRC, typically only **2–4 IBOs change substantially**. For each tracked IBO, compute the change in its per-atom IAO charge distribution between the R and P ends.
  - Bond(A–B) → lone pair(B): arrow from bond to atom B.
  - Lone pair(A) → bond(A–C): arrow from atom A to the A–C midpoint.
  - π(A=B) → σ(B–C): arrow from bond to bond.
- **Answer to "does it deduce electron movement correctly?"**: For a closed-shell step **yes, automatically**. The arrows are a direct readout of the converged wavefunction along a verified IRC, not a heuristic. What the QM **cannot** decide on its own is *which mechanism* the user meant: it computes the step that was drawn. Competing mechanisms (e.g., polar vs SET Grignard) must be drawn and compared as separate steps.
- **Teaching feature (optional)**: the user draws their expected arrows, and the system overlays "your arrows vs the wavefunction's arrows". Ketcher has no native electron-pushing arrows bound to atoms/bonds, so this needs a small custom overlay.
- **Acceptance Criteria**: For SN2, Diels–Alder and Grignard C–C formation, the auto-generated arrows match the textbook arrows.

## Phase 7 (New): Multi-Step Mechanism Assembly ("Gluing")
- **Input**: an ordered list of elementary steps, each drawn by the user: `R → I1`, `I1 → I2`, … `In → P`. Each step is marked "has TS" or "barrierless".
- **Global atom index space**: every atom that appears in *any* step (including solvent ligands, counter-ions and workup reagents) gets one permanent index. Atoms not participating in a step are **spectators**: present in the calculation, or parked far away and faded in the UI when they enter or leave (e.g., H₃O⁺ during workup).
- **Glue conditions** (checked automatically):
  1. The product end of step *k* and the reactant end of step *k+1* have the same atom set and the same connectivity.
  2. Both ends relax to the **same intermediate minimum**. That shared minimum is optimized once and used as the joint, and the frames are Kabsch-aligned on it.
  3. Energies are referenced to the separated reactants, at one level of theory.
- **Outputs**: one continuous trajectory with step boundaries marked on the timeline, a **reaction energy profile** (R, TS1, I1, TS2, … P) next to the 3D view, and orbital/arrow tracks per step.
- **Automatic intermediate discovery** (enumerating plausible intermediates from R and P alone) is research-grade and expensive. It is **out of scope**: the user supplies intermediates.

---

## Case Study: Grignard Addition (MeMgBr + H₂C=O)
Chosen to stress-test the plan: main-group metal, coordinated intermediates, solvent dependence, multiple steps.

| Step | Chemistry | Route | Notes |
|---|---|---|---|
| 0 | R–Br + Mg(0) → RMgBr | ❌ Tier D | Radical, on a metal surface. Start from the pre-formed reagent. |
| 1 | Schlenk equilibrium / aggregation | Static structures only | Show as optimized structures; real speciation is solvent-dynamic (cf. Peltzer, Eisenstein, Cascella et al., JACS 2020). |
| 2 | MeMgBr(OMe₂)₂ + H₂C=O → MeMgBr(OMe₂)(O=CH₂) + OMe₂ | Barrierless / relaxed scan | Carbonyl coordination to Mg (ligand exchange). Likely no clean TS. |
| 3 | 4-centre TS: C–C forms, Mg–C breaks, Mg–O tightens | Tier B fallback (xTB + NEB + Sella) | The key step. Polar and closed-shell, so RKS + IBO arrows work: σ(Mg–C) pair → new σ(C–C), π(C=O) → O lone pair. |
| 4 | Alkoxide + H₃O⁺ → alcohol | Barrierless / scan or cartoon | Charged and solvent-mediated. Low educational value to compute rigorously. |

- **Model system**: Me₂O instead of THF and Me instead of larger R. MeMgBr + H₂CO + 2 Me₂O ≈ 28 atoms, comfortably inside Phase 5 limits. MgCl is cheaper than MgBr if needed.
- **Tool compatibility**: RitS, OA-ReactDiff, MACE-OFF and ANI-2x **all lack Mg**, so this case *requires* the xTB/NEB fallback (or an OMol25-class model). Use ALPB(THF) implicit solvent in xTB and PCM in PySCF.
- **Lesson for the plan**: the "pure ML" pipeline covers organic textbook steps. Anything with metals depends on the physics fallback, which must be a first-class route rather than an afterthought.

---

## Human-in-the-Loop Points
The pipeline should run unattended for Tier A, but the user has defined decision points:

1. **Mechanism definition** (mandatory for multi-step): draw each elementary step, set charge and spin, and mark barrierless steps.
2. **Mapping override**: accept or fix the atom map (Ketcher mapping tool) when the Phase 1 image looks wrong.
3. **TS candidate choice**: when several distinct saddles are found, pick one (default: lowest barrier).
4. **Graph-check failure**: if the IRC lands on unexpected minima, the system shows what it found. The user either accepts it (a discovered alternative) or retries.
5. **Electron flow**: **no intervention needed**; it is derived. Optional comparison with user-drawn arrows.

---

## Open Questions / Risks
- **Hardware**: RitS wants a CUDA GPU. Is one available on the server? Otherwise the NEB fallback becomes the primary route (slower but CPU-only).
- **Turnaround**: a single step with orbitals is minutes to hours, not the ~10–30 s of current single-molecule jobs. This needs a job queue, progress reporting and a **precomputed library** of classroom mechanisms.
- **Payload size**: ~50 frames × several orbitals × mesh JSON. Consider coarser grids or level-of-detail for scrubbing, and full resolution on pause.
- **Level-of-theory consistency**: geometry (MLIP/xTB) vs orbitals (B3LYP). Document this in `scientific_choices.md` once implemented.
- **Licences**: MACE-OFF (academic licence), UMA (FAIR chemistry licence) — check before deployment.

## Suggested Order of Work
1. Phase 1 (+ H completion) → Phase 2b (xTB/NEB fallback) → Phases 3–4 with xTB. This gives a working end-to-end path for **any** element without ML dependencies.
2. Phase 5 (IBO tracking) → Phase 6 (arrows) on Tier A test set: SN2, Diels–Alder, Claisen, HCN→HNC.
3. Phase 2 primary (RitS) as an accelerator, once GPU availability is confirmed.
4. Phase 7 (multi-step), validated on the Grignard case study.

## Useful Links:
* https://github.com/chicago-aiscience/generative_chem_reactive_structures/
* https://github.com/isayevlab/RitS
* Sella (TS optimizer + IRC): https://github.com/zadorlab/sella
* Knizia & Klein, Angew. Chem. Int. Ed. 2015, 54, 5518 — electron flow from IBOs along IRCs
