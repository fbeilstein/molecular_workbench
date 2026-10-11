import os
import tempfile
import numpy as np
from .common import _unit


def _compute_vibrations(atoms):
    from ase.vibrations import Vibrations
    with tempfile.TemporaryDirectory() as tmp:
        vib = Vibrations(atoms, name=os.path.join(tmp, "vib"), delta=0.005)
        vib.run()
        data = vib.get_vibrations()
        return data.get_frequencies(), data.get_modes()


def analyse_mode(atoms, model, log, max_atoms=60):
    """Vibrational mode analysis: checks for a single imaginary frequency on reacting bonds."""
    if len(atoms) > max_atoms:
        log(f"  (skipping frequency check: {len(atoms)} atoms)")
        return None

    freqs, modes = _compute_vibrations(atoms)
    imag = [(k, f.imag) for k, f in enumerate(freqs) if f.imag > 50.0]
    report = {"n_imag": len(imag), "imag_freqs": [round(f, 1) for _, f in imag]}

    if not imag:
        log("  No imaginary frequency > 50i cm-1: this is a MINIMUM, not a TS.")
        return report

    log("  Imaginary frequencies: " + ", ".join(f"{f:.0f}i" for _, f in imag) + " cm-1")
    k0 = imag[int(np.argmax([f for _, f in imag]))][0]
    m = modes[k0]
    pos = atoms.get_positions()

    def stretch(i, j):
        return float(np.dot(_unit(pos[j] - pos[i]), m[j] - m[i]))

    all_bonds = model.r_bonds | model.p_bonds
    total = sum(abs(stretch(i, j)) for i, j in all_bonds) + 1e-12
    ch_f = [stretch(i, j) for i, j in model.formed]
    ch_b = [stretch(i, j) for i, j in model.broken]
    share = (sum(map(abs, ch_f)) + sum(map(abs, ch_b))) / total

    signs_ok = True
    if ch_f and ch_b:
        sf, sb = np.sign(ch_f), np.sign(ch_b)
        signs_ok = bool(np.all(sf == sf[0]) and np.all(sb == -sf[0]))

    report.update({"mode_share": round(share, 2), "mode_consistent": signs_ok})
    log(f"  Main imaginary mode: {share * 100:.0f}% of bond-stretch amplitude on the "
        f"reacting bonds; forming/breaking in antiphase: {'yes' if signs_ok else 'NO'}")
    return report
