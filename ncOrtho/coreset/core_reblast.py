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

# Core-set reciprocal BLAST search:
#   BLAST reference miRNA against syntenic regions in core species,
#   verify reciprocity, and build a Stockholm structural alignment.

import logging
import os
import subprocess as sp
from typing import Dict, List, Optional, Set, Tuple

try:
    from utils import check_blastdb, make_blastndb
except ModuleNotFoundError:
    from ncOrtho.utils import check_blastdb, make_blastndb


logger = logging.getLogger("ncortho")


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------
def _ranges_overlap(s1: int, e1: int, s2: int, e2: int) -> bool:
    """Return True if closed intervals [s1, e1] and [s2, e2] overlap."""
    return s1 <= e2 and s2 <= e1


def _vprint(msg: str, verbose: bool) -> None:
    """Print *msg* only when *verbose* is True."""
    if verbose:
        print(msg, flush=True)


# ---------------------------------------------------------------------------
# T-Coffee alignment
# ---------------------------------------------------------------------------
def make_alignment(
    out: str,
    mirna: str,
    cpu: int,
    core_fasta: str,
    rcoffee: str,
) -> str:
    """
    Build a ClustalW alignment with T-Coffee (optionally R-Coffee),
    then convert to Stockholm format with secondary-structure annotation.

    Returns the path to the Stockholm file.
    """
    alignment = os.path.join(out, f"{mirna}.aln")
    stockholm = os.path.join(out, f"{mirna}.sto")

    # Step 1: sequence alignment
    tc_cmd = [
        "t_coffee",
        "-quiet",
        f"-multi_core={cpu}",
    ]
    if rcoffee == "yes":
        tc_cmd.append("-mode=rcoffee")
    tc_cmd += ["-in", core_fasta, "-output=clustalw_aln", f"-outfile={alignment}"]

    sp.run(tc_cmd, capture_output=True, check=False)

    # Step 2: add structure and convert to Stockholm
    reformat_cmd = [
        "t_coffee",
        "-other_pg", "seq_reformat",
        "-in", alignment,
    ]
    if rcoffee == "yes":
        reformat_cmd += ["-action", "+add_alifold"]
    reformat_cmd += ["-output", "stockholm_aln", "-out", stockholm]

    sp.run(reformat_cmd, capture_output=True, check=False)

    return stockholm


# ---------------------------------------------------------------------------
# Reference bitscore
# ---------------------------------------------------------------------------
def maximum_blast_bitscore(
    mirna: str,
    seq: str,
    blastdb: str,
    cpu: int,
    dust: str,
) -> float:
    """
    BLAST the reference miRNA against its own genome to determine the
    maximum achievable bitscore (used as the denominator for cutoffs).

    Returns 0.0 with a warning if the sequence is not found.
    """
    cmd = [
        "blastn",
        "-num_threads", str(cpu),
        "-dust", dust,
        "-task", "megablast",
        "-db", blastdb,
        "-outfmt", "6 bitscore",
    ]
    result = sp.run(cmd, input=seq, capture_output=True, text=True, check=False)

    if result.returncode != 0 or not result.stdout.strip():
        logger.warning(
            "Reference sequence of %s not found in reference genome. "
            "Setting maximum bitscore to 0.",
            mirna,
        )
        return 0.0

    try:
        return float(result.stdout.strip().split("\n")[0].split("\t")[0])
    except (ValueError, IndexError):
        logger.warning(
            "Could not parse bitscore for %s. Setting maximum bitscore to 0.",
            mirna,
        )
        return 0.0


# ---------------------------------------------------------------------------
# Core BLAST + reciprocal verification
# ---------------------------------------------------------------------------
def _blast_synteny_regions(
    preseq: str,
    synteny_regs: str,
    cpu: int,
    dust: str,
) -> List[List[str]]:
    """
    BLAST the reference pre-miRNA against the syntenic-region database.

    Returns a list of hit rows: [sseqid, evalue, bitscore, sseq].
    """
    if not check_blastdb(synteny_regs):
        make_blastndb(synteny_regs, synteny_regs)

    cmd = [
        "blastn",
        "-num_threads", str(cpu),
        "-task", "blastn",
        "-dust", dust,
        "-db", synteny_regs,
        "-outfmt", "6 sseqid evalue bitscore sseq",
    ]
    result = sp.run(cmd, input=preseq, capture_output=True, text=True, check=False)
    if result.returncode != 0:
        raise sp.SubprocessError(
            f"BLASTn against synteny regions failed: {result.stderr.strip()}"
        )

    hits = []
    for line in result.stdout.strip().split("\n"):
        if line:
            hits.append(line.split())
    return hits


def _reblast_candidate(
    degap_seq: str,
    ref_blastdb: str,
    cpu: int,
    dust: str,
) -> List[Tuple[str, int, int]]:
    """
    Re-BLAST a candidate sequence against the reference genome.

    Returns a list of (chrom, start, end) tuples.
    """
    cmd = [
        "blastn",
        "-num_threads", str(cpu),
        "-task", "blastn",
        "-dust", dust,
        "-db", ref_blastdb,
        "-outfmt", "6 sseqid sstart send",
    ]
    result = sp.run(cmd, input=degap_seq, capture_output=True, text=True, check=False)

    reblast_hits = []
    for line in result.stdout.strip().split("\n"):
        if not line:
            continue
        fields = line.split()
        reblast_hits.append((fields[0], int(fields[1]), int(fields[2])))
    return reblast_hits


# ---------------------------------------------------------------------------
# Main entry point
# ---------------------------------------------------------------------------
def blastsearch(
    mirna: list,
    r_path: str,
    o_path: str,
    cpu: int,
    dust: str,
    verbose: bool,
    coffee: str,
) -> Optional[str]:
    """
    Perform reciprocal BLAST search for a single miRNA and construct a
    Stockholm alignment of verified core orthologs.

    Parameters
    ----------
    mirna : list
        At least 6 elements: [mirid, chrom, start, end, strand, seq, ...].
    r_path : str
        Path to the reference genome FASTA.
    o_path : str
        Base output/tmp directory.
    cpu : int
        Number of BLAST threads.
    dust : str
        ``"yes"`` or ``"no"`` — BLASTn dust filter setting.
    verbose : bool
        Print extra diagnostic output.
    coffee : str
        ``"yes"`` for r_coffee, ``"no"`` for default t_coffee.

    Returns
    -------
    str or None
        Path to the Stockholm alignment, or None on failure.
    """
    mirid, rawchrom, mstart, mend, mstrand, rawseq = mirna[:6]
    mchr = rawchrom.replace("chr", "")
    mstart = int(mstart)
    mend = int(mend)
    preseq = rawseq.replace("U", "T").replace("-", "")

    miroutdir = os.path.join(o_path, mirid)
    os.makedirs(miroutdir, exist_ok=True)
    synteny_regs = os.path.join(miroutdir, f"synteny_regions_{mirid}.fa")

    print(f"# {mirid}", flush=True)

    # If no syntenic regions were produced by the upstream pipeline,
    # train with the reference sequence alone.
    if not os.path.isfile(synteny_regs):
        logger.warning(
            "No synteny regions found for %s. Training with reference only.",
            mirid,
        )
        with open(synteny_regs, "w") as fh:
            fh.write(f">{mirid}\n{preseq}\n")
        return make_alignment(miroutdir, mirid, cpu, synteny_regs, coffee)

    # Ensure reference BLASTdb exists
    ref_basename = os.path.splitext(os.path.basename(r_path))[0]
    ref_blastdb_dir = os.path.join(o_path, "refBLASTdb")
    os.makedirs(ref_blastdb_dir, exist_ok=True)
    ref_blastdb = os.path.join(ref_blastdb_dir, ref_basename)
    if not check_blastdb(ref_blastdb):
        make_blastndb(r_path, ref_blastdb)

    max_bitscore = maximum_blast_bitscore(mirid, preseq, ref_blastdb, cpu, dust)

    # BLAST reference miRNA against syntenic regions
    ortholog_candidates = _blast_synteny_regions(preseq, synteny_regs, cpu, dust)
    _vprint(f"Number of ortholog candidates: {len(ortholog_candidates)}", verbose)

    # Reciprocal BLAST verification
    outputcol: Dict[str, str] = {}
    seen_regions: Set[str] = set()

    for hit in ortholog_candidates:
        if len(hit) < 4:
            continue
        core_region, _evalue, bitscore, sseq = hit[:4]

        if float(bitscore) <= max_bitscore * 0.5:
            continue

        degap_seq = sseq.replace("-", "")

        reblast_hits = _reblast_candidate(degap_seq, ref_blastdb, cpu, dust)
        for refchrom, refstart, refend in reblast_hits:
            if refchrom != mchr:
                _vprint(
                    f"Hit on chromosome {refchrom}, expected {mchr}",
                    verbose,
                )
                continue

            # Ensure start <= end for overlap check
            rs, re = min(refstart, refend), max(refstart, refend)
            if _ranges_overlap(rs, re, mstart, mend):
                _vprint("Reciprocity fulfilled.", verbose)
            else:
                _vprint("Reciprocity unfulfilled.", verbose)
                continue

            # Only keep the best hit per syntenic region
            if core_region not in seen_regions:
                outputcol[core_region] = degap_seq
                seen_regions.add(core_region)

    # Write core orthologs FASTA
    corefile = os.path.join(miroutdir, f"core_orthologs_{mirid}.fa")
    with open(corefile, "w") as fh:
        fh.write(f">reference\n{preseq}\n")
        if outputcol:
            for region, sequence in outputcol.items():
                fh.write(f">{region}\n{sequence}\n")
        else:
            logger.warning(
                "No core orthologs found for %s. Training with reference only.",
                mirid,
            )

    return make_alignment(miroutdir, mirid, cpu, corefile, coffee)