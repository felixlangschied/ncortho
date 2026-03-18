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

# Reverse-BLAST verification of CM/pHMM candidate orthologs.

import logging
import os
import subprocess as sp
from typing import Optional

try:
    from blastparser import BlastParser
except ModuleNotFoundError:
    from ncOrtho.blastparser import BlastParser


logger = logging.getLogger("ncortho")


def perform_reblast(
    sequence: str,
    refblast: str,
    cpu: int,
    outdir: str,
    candidate: str,
    mirna_data,
    dust: str,
    msl: float,
    cleanup: bool,
    check_coorth_ref: bool,
) -> str:
    """
    Run a reverse BLASTn of *sequence* against the reference database
    and check whether the best hit overlaps the original miRNA.

    Parameters
    ----------
    sequence : str
        Nucleotide sequence of the candidate ortholog.
    refblast : str
        Path to the reference BLASTdb.
    cpu : int
        Number of threads for BLASTn.
    outdir : str
        Directory for intermediate files.
    candidate : str
        Candidate identifier (used for filenames).
    mirna_data : RNA object
        Reference miRNA with positional attributes.
    dust : str
        ``"yes"`` or ``"no"`` — BLASTn dust filter setting.
    msl : float
        Minimum sequence length ratio.
    cleanup : bool
        If False, keep the raw BLAST output file.
    check_coorth_ref : bool
        If True and the best hit does not overlap the reference, check
        whether it is a co-ortholog.

    Returns
    -------
    str
        The candidate *sequence* if accepted, otherwise an empty string.
    """
    # --- Run BLASTn --------------------------------------------------------
    blast_cmd = [
        "blastn",
        "-task", "blastn",
        "-db", refblast,
        "-num_threads", str(cpu),
        "-dust", dust,
        "-outfmt", "6 sseqid sstart send sstrand bitscore",
    ]
    result = sp.run(
        blast_cmd,
        input=sequence,
        capture_output=True,
        text=True,
        check=False,
    )
    if result.returncode != 0:
        logger.warning("re-BLAST failed for %s: %s", candidate, result.stderr.strip())
        return ""

    raw_output = result.stdout

    # Optionally persist raw output for debugging
    if not cleanup:
        blast_output_path = os.path.join(outdir, f"reBLAST_{candidate}.out")
        with open(blast_output_path, "w") as fh:
            fh.write(raw_output)

    # --- Parse results -----------------------------------------------------
    blast_lines = [
        line.split()
        for line in raw_output.strip().split("\n")
        if line.strip()
    ]
    if not blast_lines:
        logger.info("No re-BLAST hits for %s", candidate)
        return ""

    # --- Evaluate reciprocal best hit --------------------------------------
    bp = BlastParser(mirna_data, blast_lines, msl)

    if bp.evaluate_besthit():
        logger.info("Accepted %s: best hit overlaps reference miRNA", candidate)
        return sequence

    if check_coorth_ref:
        logger.info(
            "Best hit for %s differs from reference; checking co-ortholog status",
            candidate,
        )
        if bp.check_coortholog_ref(sequence, outdir):
            return sequence

    logger.info("Rejected %s: best hit does not overlap reference miRNA", candidate)
    return ""