# Transition State Search Strategies, Reaction Coordinate Scanning, and Leaving Groups

## 1. Executive Summary
This document provides the foundational theory, algorithms, failure modes, and architectural fixes for Transition State (TS) optimization in the Reaction Workbench.
Special focus is given to **unimolecular elimination / bond-cleaving reactions** (such as tetrahedral intermediate decomposition in nucleophilic acyl substitution: $\text{R-C(O}^-)(\text{Cl})\text{-R}' \to \text{R-CO-R}' + \text{Cl}^-$).

---

## 2. The Nature of the Problem: Why Generic TS Guesses Fail on Leaving Groups

### 2.1 The Hammond Postulate and Early Transition States
In bimolecular additions ($\text{A} + \text{B} \to [\text{A}\cdots\text{B}]^\ddagger$), forming bonds typically elongate by $+0.5\text{--}0.8\text{ \AA}$ relative to covalent $d_0$.
Similarly, in homolytic cleavages or difficult eliminations ($\text{C}-\text{C}$ or $\text{C}-\text{H}$ cleavage), breaking bonds elongate by $+0.35\text{--}0.50\text{ \AA}$ before the saddle point is reached.

However, for **exceptional leaving groups** (halides like $\text{Cl}^-$, $\text{Br}^-$, sulfonates like $\text{OTs}^-$, $\text{OTf}^-$) assisted by an adjacent electron donor (such as an oxyanion $\text{O}^-$ reforming a carbonyl $\text{C}=\text{O}$):
1. The driving force into products is large and exergonic.
2. The transition state is **extremely early** (Hammond's postulate).
3. The potential energy well of the tetrahedral intermediate is shallow ($0.5\text{--}4.5\text{ kcal/mol}$), and the saddle point occurs at a very specific distance along the $\text{C}-\text{X}$ coordinate.

### 2.2 The Root Failure Modes in the Legacy Codebase
1. **Static Target Heuristic ($d_0 + 0.35\text{ \AA}$):**
   * Legacy routines in `lego_builder.py`, `manual_valley.py`, `pose_sampler.py`, and `rits_predictor.py` computed breaking bond targets as $d_{\text{target}} = d_0(i, j) + 0.35\text{ \AA}$.
   * For $\text{C}-\text{Cl}$ ($d_0 \approx 1.78\text{ \AA}$), this forced the initial guess to $2.13\text{ \AA}$.
   * At $2.13\text{ \AA}$, the $\text{C}-\text{Cl}$ bond is already on the steep downhill product slope (several $\text{kcal/mol}$ below the barrier).
2. **Unconstrained Sella Divergence:**
   * Starting Sella from $d > 2.08\text{ \AA}$ placed it in an outward gradient of $> 1.1\text{ eV/\AA}$.
   * Following the eigenvector outward caused Sella to drift to $2.42\text{ \AA}$ (toward separated products) instead of climbing backward to the saddle point.
3. **Premature IRC Reverse Abort (`irc_runner.py`):**
   * Integrating backwards (in reverse) from $2.08\text{ \AA}$ toward the intermediate at $1.78\text{ \AA}$ required climbing uphill back to the saddle point.
   * A check `if e > min(raw) + 0.04 eV: break` aborted the reverse integration after only $0.01\text{ \AA}$, mistaking the uphill climb toward the saddle for leaving a basin!
   * As a result, the "reactant" in Step 2 was falsely pinned at $2.07\text{ \AA}$, giving a forward barrier of $0.00\text{ kcal/mol}$.

---

## 3. Systematic Audit of Step 3 Strategies

| Strategy | Description | Applicability to Bond Breaking / Elimination | Status & Required Fix |
| :--- | :--- | :--- | :--- |
| **`drive` (Constrained 1D Scan)** | Scans the reaction coordinate incrementally with orthogonal relaxation, finds the energy maximum, and refines. | **100% Ideal (Primary Strategy)** | **Fixed:** Start from bound minimum ($d_0 - 0.05\text{ \AA}$), scan in $0.02\text{--}0.03\text{ \AA}$ steps across the barrier window, and guard against Sella drift. |
| **`lego` (Dual-Endpoint IDPP)** | Interpolates reactant and product conformers via IDPP, then polishes with Sella. | **Applicable** | **Fixed:** Do not hardcode $+0.35\text{ \AA}$ for breaking bonds; perform a multi-point IDPP ridge search. |
| **`manual_valley` (Quantum Valley)** | Relaxes crude manual placements with harmonic restraints. | **Applicable (with manual positioning)** | **Fixed:** Preserve the user's manual distance rather than overwriting with $d_0 + 0.35\text{ \AA}$. |
| **`neb` (Climbing Image NEB)** | Nudged Elastic Band connecting R and P. | **Applicable** | **Fix:** Increase image density in the near-barrier region so narrow barriers ($< 0.1\text{ \AA}$) are not stepped over. |
| **`orca_dft`** | High-level DFT (r2SCAN-3c) transition state optimization and validation. | **Applicable** | Requires a well-conditioned starting guess from `drive` or `lego`. |
| **`genetic_swarm`** | SE(3) rigid-body genetic algorithm. | **Not Applicable** | Designed exclusively for bimolecular multi-fragment docking ($\ge 2$ reactants). No effect on single-molecule elimination. |
| **`auto_dock`** | 32-pose geometric approach sampler. | **Not Applicable** | Exclusively for bimolecular approach orientations. |
| **`dimer`** | Min-mode following without product endpoint. | **Marginal** | Often overshoots shallow barriers on semiempirical surfaces. |
| **`rits`** | Generative Flow Matching (ML). | **Conditional** | Requires external neural checkpoint (`rits.ckpt`). |

---

## 4. Architectural Solutions Implemented

### 4.1 Fine-Grid Constrained Coordinate Drive (`run_drive_sella_ts`)
1. For bond-breaking steps, the coordinate drive begins at $d_{\text{start}} = \min(d_{\text{current}}, d_0 - 0.05\text{ \AA})$.
2. The distance is driven outward in fine increments ($\Delta d \approx 0.025\text{ \AA}$) up to $d_0 + 0.50\text{ \AA}$.
3. At each step, all orthogonal degrees of freedom are fully relaxed via LBFGS while restraining the reacting bond.
4. The exact maximum along the energy profile $E(d)$ is identified as the saddle candidate.
5. If Sella is used to polish the saddle point, a tight trust-zone guard ($\text{max\_drift} \le 0.05\text{ \AA}$) prevents runaway divergence into the product basin.

### 4.2 Interactive 3D Reacting Bond Distance Control (UI)
When reacting bonds are selected in the 3D TS Editor:
* The user is provided with a real-time **Bond Distance Slider and Numeric Input (in Å)**.
* Adjusting the slider dynamically repositions the leaving group atom along the internuclear axis in the 3D viewport.
* This allows the user to inspect and manually set the exact distance (e.g. $1.85\text{ \AA}$) before running optimization or sending the geometry to IRC.

### 4.3 Robust Reverse IRC Descent
In `irc_runner.py`:
* When integrating backwards into a bound unimolecular reactant, the descent vector contracts the reacting bond toward the covalent minimum $d_0$.
* The integration terminates only when the bond reaches covalent equilibrium ($d \le d_0$) or gradient norm drops below $f_{\text{max}}$, preventing premature abort on minor slope fluctuations.

### 4.4 Universal Physics-Based 1D Coordinate Drive Scaling
In `reaction_model.py`, `strategies.py`, and `ts_editor_cards.js`:
* **Zero Hardcoding**: Scan ranges and grid points are never hardcoded. They are derived directly from the fundamental covalent radii ($d_0 = r_{\text{cov}}(A) + r_{\text{cov}}(B)$) and current 3D interatomic distance ($d_{\text{cur}}$) of the reacting pair:
  * **Bond-Cleaving Reactions** ($(i, j) \in \text{model.broken}$):
    $$d_{\text{start}} = \begin{cases} d_{\text{cur}}, & \text{if } d_{\text{cur}} < 1.35 \times d_0 \\ d_0, & \text{otherwise} \end{cases}$$
    $$d_{\text{end}} = \max\left(d_{\text{start}} + 0.42 \times d_0, \; 1.55 \times d_0\right)$$
    - For C-H cleavage ($d_0 = 1.07\text{ \AA}$): scans $1.07\text{ \AA} \to 1.66\text{ \AA}$ (16 points).
    - For C-C cleavage ($d_0 = 1.52\text{ \AA}$): scans $1.52\text{ \AA} \to 2.36\text{ \AA}$ (22 points).
    - For C-Cl cleavage ($d_0 = 1.78\text{ \AA}$): scans $1.78\text{ \AA} \to 2.76\text{ \AA}$ (25 points).
    - For elongated/zwitterionic wells (e.g. C-Cl at $2.06\text{ \AA}$): dynamically scans $2.06\text{ \AA} \to 2.81\text{ \AA}$ (20 points).
  * **Bond-Forming Reactions** ($(i, j) \in \text{model.formed}$):
    $$d_{\text{start}} = \min\left(\max(d_{\text{cur}}, d_0 + 1.0), \; d_0 + 2.5\right)$$
    $$d_{\text{end}} = d_0$$
  * **Grid Resolution**:
    $$N_{\text{points}} = \text{clamp}\left(\text{round}\left(\frac{|d_{\text{end}} - d_{\text{start}}|}{0.04}\right) + 1, \; 15, \; 40\right)$$
* **Dynamic Prefill & Customization**:
  * The TS Editor UI dynamically computes and displays these physical parameters for the currently active reaction and live 3D coordinates via `/api/step/step_3_ts/custom` (`compute_scan_defaults`).
  * The `↺ Prefill Defaults` button dynamically recalculates from the live 3D coordinates in the viewport.
  * If the user enters custom values, they override the defaults. If inputs are left blank, the backend automatically evaluates the physical defaults.
  * Rigid geometric constraint via ASE `FixBondLength` ensures that each point is maintained to exact numerical precision during orthogonal relaxation.


