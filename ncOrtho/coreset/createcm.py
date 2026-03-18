"""
ncOrtho - Targeted ortholog search for miRNAs
Copyright (C) 2021 Felix Langschied

ncOrtho is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

ncOrtho is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with ncOrtho.  If not, see <http://www.gnu.org/licenses/>.
"""

# Construct and calibrate covariance models (CMs) or profile HMMs
# for ncRNA core-set alignments in Stockholm format.

import logging
import os
import subprocess as sp
from typing import Optional


logger = logging.getLogger("ncortho")


# ---------------------------------------------------------------------------
# CM construction and calibration
# ---------------------------------------------------------------------------
def _run_cmd(cmd: list, description: str) -> sp.CompletedProcess:
    """Run a command, raising on failure with a clear message."""
    result = sp.run(cmd, capture_output=True, text=True, check=False)
    if result.returncode != 0:
        raise RuntimeError(
            f"{description} failed (exit {result.returncode}): "
            f"{result.stderr.strip()}"
        )
    return result


def create_cm(
    alignment: str,
    outpath: str,
    name: str,
    cpu: int,
) -> str:
    """
    Build and calibrate a covariance model from a Stockholm alignment.

    Parameters
    ----------
    alignment : str
        Path to the input alignment (``.sto``).
    outpath : str
        Directory for the output model.
    name : str
        Model name (used for ``-n`` and the output filename).
    cpu : int
        Number of CPU cores for calibration.

    Returns
    -------
    str
        Path to the calibrated ``.cm`` file.
    """
    model = os.path.join(outpath, f"{name}.cm")

    build_cmd = [
        "cmbuild",
        "-F",
        "-n", name,
        model,
        alignment,
    ]
    _run_cmd(build_cmd, f"cmbuild for {name}")

    calibrate_cmd = [
        "cmcalibrate",
        "--cpu", str(cpu),
        model,
    ]
    _run_cmd(calibrate_cmd, f"cmcalibrate for {name}")

    return model


def create_phmm(
    alignment: str,
    outpath: str,
    name: str,
    cpu: int,
) -> str:
    """
    Build a profile HMM from a Stockholm alignment.

    Parameters
    ----------
    alignment : str
        Path to the input alignment (``.sto``).
    outpath : str
        Directory for the output model.
    name : str
        Model name.
    cpu : int
        Number of CPU cores (currently unused by hmmbuild but kept
        for API consistency with ``create_cm``).

    Returns
    -------
    str
        Path to the ``.phmm`` file.
    """
    model = os.path.join(outpath, f"{name}.phmm")

    cmd = [
        "hmmbuild",
        "-n", name,
        "--informat", "stockholm",
        "--dna",
        model,
        alignment,
    ]
    _run_cmd(cmd, f"hmmbuild for {name}")

    return model