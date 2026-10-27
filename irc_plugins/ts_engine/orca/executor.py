"""
ORCA Executor Module
Handles subprocess execution, environment configuration, and basic output extraction.
"""

import os
import subprocess
import shutil
import tempfile
import re
from typing import Any, Dict, List, Optional

ORCA_BIN = "/data/orca/orca"
os.environ["OMPI_MCA_rmaps_base_oversubscribe"] = "1"
os.environ["OMPI_MCA_btl_vader_single_copy_mechanism"] = "none"


def run_orca_calc(
    symbols: List[str],
    positions: Any,
    charge: int,
    mult: int,
    solvent: str,
    extra_keywords: str,
    block_settings: str = "",
    log: Optional[Any] = None,
    tmpdir: Optional[str] = None,
    name: str = "calc",
    nprocs: int = 4
) -> Dict[str, Any]:
    """Run an ORCA calculation and return parsed output dictionary."""
    cleanup = False
    if tmpdir is None:
        tmpdir = tempfile.mkdtemp(prefix=f"orca_{name}_")
        cleanup = True

    inp_file = os.path.join(tmpdir, f"{name}.inp")
    out_file = os.path.join(tmpdir, f"{name}.out")
    solv_spec = f"CPCM({solvent.capitalize()})" if solvent and solvent.lower() != "gas" else ""

    with open(inp_file, "w") as f:
        f.write(f"! r2SCAN-3c {solv_spec} {extra_keywords}\n")
        f.write(f"%pal nprocs {nprocs} end\n")
        if block_settings:
            f.write(f"{block_settings}\n")
        f.write(f"* xyz {charge} {mult}\n")
        for sym, pos in zip(symbols, positions):
            f.write(f"  {sym:<2} {pos[0]:12.6f} {pos[1]:12.6f} {pos[2]:12.6f}\n")
        f.write("*\n")

    if log:
        log(f"  [ORCA 6.1] Running {extra_keywords} ({name})...")

    with open(out_file, "w") as out_fh:
        proc = subprocess.run([ORCA_BIN, inp_file], stdout=out_fh, stderr=subprocess.STDOUT, cwd=tmpdir)

    out_txt = ""
    if os.path.exists(out_file):
        with open(out_file, "r") as f:
            out_txt = f.read()

    res = {
        "success": proc.returncode == 0,
        "out_file": out_file,
        "tmpdir": tmpdir,
        "output": out_txt,
        "energy_eh": None,
        "gibbs_eh": None,
        "zpe_kcal": None,
        "imaginary_freqs": []
    }

    for line in out_txt.splitlines():
        if "FINAL SINGLE POINT ENERGY" in line:
            try:
                res["energy_eh"] = float(line.split()[-1])
            except ValueError:
                pass
        elif "Final Gibbs free energy" in line:
            m = re.search(r"Final Gibbs free energy\s+\.\.\.\s+([-+]?\d*\.?\d+)", line)
            if m:
                res["gibbs_eh"] = float(m.group(1))
        elif "Zero point energy" in line and "kcal/mol" in line:
            m = re.search(r"(\d*\.?\d+)\s+kcal/mol", line)
            if m:
                res["zpe_kcal"] = float(m.group(1))
        elif "***imaginary mode***" in line:
            m = re.search(r"([-+]?\d*\.?\d+)\s+cm\*\*-1", line)
            if m:
                res["imaginary_freqs"].append(float(m.group(1)))

    if cleanup:
        shutil.rmtree(tmpdir, ignore_errors=True)
    return res
