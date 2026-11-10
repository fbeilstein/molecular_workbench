# Reaction Free Energy and Solution Thermodynamics

This document outlines the theoretical foundation, mathematical formulation, and software architecture for computing **Gibbs Free Energy ($\Delta G^\circ_{solv}$)** profiles along Intrinsic Reaction Coordinate (IRC) trajectories and multistep reaction cascades.

---

## 1. Physical Foundation: Electronic Energy vs. Free Energy

When electronic structure engines (such as xTB or ORCA DFT) compute single points or downhill gradients along an IRC trajectory, they evaluate the **Born-Oppenheimer Electronic Potential Energy ($E_{elec}$)**.

In experimental chemistry, reactions occur at standard state in liquid solution ($T = 298.15\text{ K}$, $C^\circ = 1\text{ M}$, $p^\circ = 1\text{ atm}$), governed by the **Gibbs Free Energy**:

$$\Delta G^\circ_{solv} = \Delta E_{elec} + \Delta G_{solv} + \Delta G_{vib} - T\Delta S_{trans/rot}^\circ$$

### The Bimolecular Association Entropy Paradox on Pure Potential Energy Surfaces
On a pure electronic potential energy surface $E_{elec}$:
- Charged or polar reactants (e.g. $\text{CH}_3\text{COCl} + \text{HCOO}^-$) experience strong attractive electrostatic and dispersion interactions upon approach.
- As a consequence, $E_{elec}(\text{TS})$ is frequently **lower** than the sum of separated reactants $E_{elec}(A) + E_{elec}(B)$, falsely indicating a "negative activation barrier" or barrierless process.
- In reality, condensing two independent free molecules in solution into a single transition state complex $[AB]^\ddagger$ incurs a large loss of 3 translational and 3 rotational degrees of freedom (converted into weak intermolecular vibrations/librations).

According to the classical thermodynamic models of Whiting & Pacey (1974) and Amzel & Eckstein (1979) for standard state $1\text{ M}$ solution:
$$\Delta G_{assoc}^\circ \approx +10 \dots +14\text{ kcal/mol} \quad (\text{default standard value: } 13.76\text{ kcal/mol})$$

This entropic penalty converts the apparently negative electronic barrier into the experimentally observed positive activation free energy $\Delta G^\ddagger > 0$.

---

## 2. Molecularity & Reaction Types

The thermodynamic engine (`irc_plugins/ts_engine/free_energy_engine.py`) automatically determines the reaction molecularity:

### Case 1: Association ($A + B \longrightarrow C$, e.g. Nucleophilic Addition)
- **Reactants ($s \to -\infty$):** $N_R = 2$ independent molecules.
- **Transition State ($s = 0$):** $N_{TS} = 1$ bound activated complex.
- **Product ($s \to +\infty$):** $N_P = 1$ bound covalent intermediate (e.g. tetrahedral adduct).
- **Thermodynamic Rule:** Both the TS and the product adduct are single bound entities. Both retain the exact same association entropy penalty ($+13.76\text{ kcal/mol}$).
- **Continuity:** Because both TS and the adduct share the same molecularity ($N=1$), **no energy cliff exists between TS and the intermediate**.

### Case 2: Dissociation ($C \longrightarrow A + B$, e.g. Leaving Group Elimination)
- **Reactant:** $N_R = 1$ bound molecule.
- **Transition State:** $N_{TS} = 1$ bound activated complex.
- **Products:** $N_P = 2$ independent departing fragments (e.g. ester + $\text{Cl}^-$).
- **Thermodynamic Rule:** Upon asymptotic separation ($d \to \infty$), the departing fragments gain translational and rotational entropy:
  $$\Delta G_{dissoc}^\circ = -(N_P - 1) \times 13.76\text{ kcal/mol}$$
- **Continuity:** As the breaking bond stretches, the free energy curve smoothly drops towards the asymptotic dissociated limit.

### Case 3: Multi-Fragment Extrusion ($A \longrightarrow B + C + D$)
- **Reactant:** $N_R = 1$ molecule.
- **Products:** $N_P = 3$ independent fragments (e.g. decarboxylative elimination).
- **Thermodynamic Rule:** Gains two quanta of dissociation entropy:
  $$\Delta G_{dissoc}^\circ = -(3 - 1) \times 13.76 = -27.52\text{ kcal/mol}$$

### Case 4: Spectators and Catalysts ($A + B \longrightarrow C + A$)
- **Coordinating Catalyst / Reactive Spectator:**
  - Enters the TS complex $[A\cdot B]^\ddagger$ ($N_{TS} = 1$) paying the association entropy penalty.
  - Leaves the product state as an independent molecule ($N_P = 2$), fully regaining its entropy.
  - In net reaction free energy $\Delta G_{rxn}^\circ$, the spectator's energy cancels out identically via Hess's Law.
- **Non-coordinating Spectator (Passive counter-ion):**
  - Stays outside the coordination sphere throughout the entire reaction; its inter-fragment distance remains $> d_{contact}$. It is automatically excluded from the activation penalty.

---

## 3. Mathematical Formulation & Continuity Guarantee

To ensure the energy curve $G(s)$ is strictly $C^1$-continuous (differentiable, without any cliffs or derivative spikes), fragment separation is governed by an analytic switching function $\lambda(d)$:

$$\lambda(d) = \begin{cases}
0, & d \le d_{contact} \quad (\text{bound complex}) \\
\sin^2\left(\frac{\pi}{2} \cdot \frac{d - d_{contact}}{d_\infty - d_{contact}}\right), & d_{contact} < d < d_\infty \\
1, & d \ge d_\infty \quad (\text{isolated free particles})
\end{cases}$$

Where:
- $d_{contact} = d_0 + 0.60\text{ \AA}$ (covalent bond contact threshold)
- $d_\infty = \max(5.40\text{ \AA}, d_{contact} + 3.20\text{ \AA})$ (asymptotic separation limit)

### Per-Frame Free Energy Equations
For any frame $k$ along the reaction coordinate:

1. **Reverse branch (Reactant approach, $k \le ts\_index$):**
   $$\Delta G(k) = (E_{raw}(k) - E_{raw}(0)) \times 23.0605 + (1 - \lambda_R(k)) \times (N_R - 1) \times \Delta G_{assoc}^\circ$$

2. **Forward branch (Product departure, $k \ge ts\_index$):**
   $$\Delta G(k) = (E_{raw}(k) - E_{raw}(0)) \times 23.0605 + (N_R - 1) \times \Delta G_{assoc}^\circ - \lambda_P(k) \times (N_P - 1) \times \Delta G_{assoc}^\circ$$

### Boundary Verification:
- **At $k = 0$ ($d \ge d_\infty \implies \lambda_R = 1$):**
  $\Delta G(0) = 0.0\text{ kcal/mol}$ (reference baseline).
- **At $k = ts\_index$ ($\lambda_R = 0, \lambda_P = 0$):**
  $$\Delta G(TS) = (E_{raw}(TS) - E_{raw}(0)) \times 23.0605 + (N_R - 1) \times \Delta G_{assoc}^\circ$$
  Both left and right equations yield the **exact same value and zero derivative**, guaranteeing zero cliff.
- **At $k = ts\_index + 1$ (bound intermediate):**
  $\lambda_P = 0 \implies$ follows the pure quantum mechanical potential energy surface anchored smoothly to TS.
- **At $k = N - 1$ ($d \ge d_\infty \implies \lambda_P = 1$):**
  $$\Delta G_{rxn}^\circ = (E_{raw}(N-1) - E_{raw}(0)) \times 23.0605 + (N_R - N_P) \times \Delta G_{assoc}^\circ$$
  Accurately matches standard-state Hess's Law reaction free energy.

---

## 4. Software Architecture

```
irc_plugins/
├── ts_engine/
│   ├── irc_runner.py          # Geometric IRC downhill path integration (PES only)
│   ├── free_energy_engine.py  # Rigorous thermodynamic ΔG calculation & continuity
│   └── energy_balancer.py     # Multistep cascade stoichiometry & Hess's Law balancing
├── step_4_irc.py              # Step 4 plugin: calls irc_runner + free_energy_engine
└── step_5_cascade.py          # Step 5 plugin: stitches multistep cascade with balanced ΔG
```

### Separation of Concerns
1. `irc_runner.py` is strictly a trajectory integrator. It does not contain any hardcoded magic numbers or heuristic spline drawings.
2. `free_energy_engine.py` receives the raw IRC trajectory, inspects chemical connectivity and stoichiometry, and produces the continuous Free Energy profile.
3. Unit tests in `tests/test_free_energy_continuity.py` enforce regression prevention on continuity across all reaction types.

