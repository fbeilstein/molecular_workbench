# Scientific Choices & Methodology

## 1. Geometry Optimization
Before performing the heavy ab initio PySCF calculations, the molecular geometry must be optimized to 3D. We provide two options via `mol_prep.py`:
- **GFN2-xTB (Default)**: A fast, semi-empirical method. It is extremely robust for organic and organometallic compounds and yields geometries very close to high-level DFT.
- **PySCF Geometry Optimization**: Using the `geometric` optimizer with the chosen level of theory (e.g., B3LYP). This offers ab-initio rigorous optimization when xTB approximations are insufficient.

## 2. Ab Initio Level of Theory
The pipeline uses Restricted Hartree-Fock (RHF) / Kohn-Sham Density Functional Theory (DFT) with the **B3LYP** functional and the **6-31G*** basis set.
- **Why B3LYP/6-31G*?** It is the classic "workhorse" of computational organic chemistry. It strikes the perfect balance between accuracy and computational speed, allowing the application to process queries in real-time (~10-30 seconds) on the web backend, while still providing polarization functions (`*`) on heavy atoms for accurate angular momentum descriptions (crucial for $\pi$ systems and proper lone pair shapes).

## 3. Antibond Projection (The $\sigma^*$ / $\pi^*$ Problem)
In an extended basis set like 6-31G*, the lowest-energy unoccupied orbitals (like the canonical LUMO) are often diffuse, Rydberg-like states, which act like giant bubbles waiting to capture free electrons. 
- **The Problem:** In introductory chemistry, the LUMO is taught as being the valence $\sigma^*$ or $\pi^*$ antibond. However, in modern basis sets, the true valence antibond is often pushed high up into the virtual continuum, completely masked by these diffuse states. Localizing the entire virtual space creates mathematical "donut" artifacts that look nothing like textbook antibonds.
- **The Solution:** We implemented an **Exact Antisymmetrization Projection** in `valence.py`. We take the beautifully localized $\sigma$ and $\pi$ bonding IBOs, mathematically subtract the atomic contributions ($C_A - C_B$), and project this strictly back into the virtual (empty) subspace. This forces the creation of textbook-perfect, strictly orthogonal valence antibonding states ($\sigma^*$ and $\pi^*$) that cleanly mirror their parent bonds, exactly as chemistry students expect.

## 4. Intrinsic Bond Orbitals (IBO) vs. Foster-Boys
We use Knizia's Intrinsic Bond Orbitals (IBO) for occupied orbital localization instead of Foster-Boys or Pipek-Mezey.
- **Why?** IBO perfectly preserves $\sigma$-$\pi$ separation. Foster-Boys minimizes the spatial extent of orbitals, which mathematically forces multiple bonds (like a $C=O$ double bond) to mix their $\sigma$ and $\pi$ components into two equivalent "banana bonds". While mathematically valid, banana bonds confuse undergraduate students. IBO correctly recovers the standard textbook picture of one tight $\sigma$ bond and one orthogonal $\pi$ bond.

## 5. The "Rabbit Ear" Lone Pair Myth
A known quirk of IBO/Pipek-Mezey localization (and modern quantum chemistry in general) is that the lone pairs on molecules like $H_2O$ and $CH_2O$ do **not** form the two equivalent $sp^3$/$sp^2$ "rabbit ears" taught in VSEPR theory. Instead, they naturally separate into:
1. A deep, $s$-rich spherical lone pair (pointing along the bond axis).
2. A higher-energy, pure $p$-character lone pair (perpendicular to the axis).
- **Scientific Choice:** We explicitly chose to preserve this physically accurate representation rather than artificially hybridizing them into rabbit ears. Photoelectron Spectroscopy (PES) definitively proves the rabbit ears don't physically exist (water has two distinct ionization energies for its lone pairs, not one). We render the true $s$- and $p$-lone pairs and allow this to be a "teachable moment" for educators to show the limits of VSEPR theory.

## 6. Population-Based Semantic Classification
We classify localized orbitals dynamically using Mulliken populations and a basis-set agnostic nodal-plane test in `classification.py`:
- **Bonds vs Lone Pairs:** If an orbital has $>10\%$ electron population on a second atom, it is a bond. This $10\%$ threshold is critical; highly polarized bonds (like $C=O$) have lone pairs with a tiny 4-5% "tail" on the Carbon atom. A cutoff lower than 10% would falsely flag these lone pairs as bonds.
- **$\sigma$ vs $\pi$:** To rigorously separate $\sigma$ and $\pi$ bonds regardless of atomic basis functions (e.g., $d$ or $f$ orbital involvement), we sample the MO on a 12-point circle perpendicular to the bond axis at its midpoint. A $\sigma$ bond is cylindrically symmetric (all points share the same sign), while a $\pi$ bond has a nodal plane passing through the axis (points exhibit a half-positive, half-negative distribution).

## 7. Extended Conjugation & Canonical Anti-Symmetry
When extracting the canonical $\pi$ system, the pipeline extracts the *entire* continuous conjugated network. 
- **Physical Rationale:** Quantum mechanics does not isolate an aromatic ring from adjacent $\pi$ bonds. Electrons delocalize over the complete conjugated topology.
- **The Anti-Symmetry Check:** To robustly identify true $\pi$ orbitals across non-planar 3D geometries, we sample points slightly above and below the calculated least-squares plane of the conjugated fragment. True $\pi$ orbitals are strictly anti-symmetric across this plane. By demanding both opposite signs and symmetric magnitudes (`abs(u - d) / max_val < 0.45`), we flawlessly filter out arbitrary $\sigma$ and lone-pair orbitals that coincidentally pass through the nuclear plane.

## 8. Hypervalency & The "Triple Bond" Anomaly
In highly polarized hypervalent oxides (e.g., the terminal oxygens in $ClO_4^-$ or $SO_4^{2-}$), the two perpendicular oxygen lone pairs are degenerate and both donate significantly (often $>10\%$) into the empty $d$-orbitals (or $\sigma^*$) of the central atom. 
- **The Problem:** A strict population threshold causes the algorithm to classify *both* of these polarized donations as $\pi$ bonds, leading to rendering impossible "triple bonds" to terminal oxygens.
- **The Heuristic Solution:** We enforce a strict chemical rule in `valence.py`: When a $\pi$ bond is detected between a hypervalent central atom (S, P, Cl, Br, I, N) and a terminal oxygen, we lock it in. If a second, highly-polarized $\pi$-donation is detected on that exact same oxygen pair, it is gracefully intercepted and demoted back to a localized oxygen lone pair. This robustly forces the output to match textbook Lewis structures (exactly one double bond per hypervalent oxo group) while preserving the underlying rigor of the IBO localization.
