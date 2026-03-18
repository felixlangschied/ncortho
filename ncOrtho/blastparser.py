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

# Reverse-BLAST result evaluation:
#   - BlastParser: verify that the best re-BLAST hit overlaps the
#     original reference miRNA (reciprocal best-hit criterion).
#   - ReBlastParser: detect co-orthologs by comparing sequence distances
#     among multiple accepted candidates.

import logging
import os
import subprocess as sp
from typing import Dict, List, Optional, Tuple

from Bio import AlignIO
from Bio.Phylo.TreeConstruction import DistanceCalculator


logger = logging.getLogger("ncortho")


# ---------------------------------------------------------------------------
# Distance-matrix helper
# ---------------------------------------------------------------------------
def _safe_remove(path: str) -> None:
    """Remove *path* if it exists; silently ignore otherwise."""
    try:
        os.remove(path)
    except FileNotFoundError:
        pass


def calculate_distance_matrix(aln_path: str):
    """
    Align sequences with t_coffee and compute an identity-based distance
    matrix.

    The input FASTA at *aln_path* is overwritten by the alignment and
    then removed (along with any ``.dnd`` tree file t_coffee creates).

    Returns
    -------
    Bio.Phylo.TreeConstruction.DistanceMatrix
    """
    env = os.environ.copy()
    env["MAX_N_PID_4_TCOFFEE"] = "4194304"

    aln_dir = os.path.dirname(os.path.abspath(aln_path))
    aln_basename = os.path.basename(aln_path)
    dnd_name = os.path.splitext(aln_basename)[0] + ".dnd"
    # t_coffee writes the .dnd file next to the input by default
    dnd_path = os.path.join(aln_dir, dnd_name)
    # It can also drop it into the cwd — cover both.
    dnd_path_cwd = os.path.join(os.getcwd(), dnd_name)

    # t_coffee's wrapper script uses bash-specific syntax (e.g. [[).
    # Running it without a shell (the default for list-based sp.run)
    # can cause "[[: not found" errors on systems where /bin/sh is dash.
    # We therefore invoke t_coffee through a shell string, matching the
    # original working invocation, but with proper quoting via shlex.
    #
    # The = syntax (-type=dna, -output=fasta_aln) must be preserved as
    # single tokens — t_coffee does not accept "-type dna" as two
    # separate arguments.
    import shlex
    cmd_str = (
        "t_coffee {infile} -no_warning -quiet -type=dna "
        "-output=fasta_aln -outfile={outfile}"
    ).format(
        infile=shlex.quote(aln_path),
        outfile=shlex.quote(aln_path),
    )

    try:
        sp.run(
            cmd_str, shell=True, env=env, check=True,
            capture_output=True, text=True,
        )
        with open(aln_path) as fh:
            aln = AlignIO.read(fh, "fasta")
        calculator = DistanceCalculator("identity")
        dm = calculator.get_distance(aln)
        return dm
    except sp.CalledProcessError as exc:
        raise RuntimeError(
            f"t_coffee alignment failed (exit {exc.returncode}): "
            f"{exc.stderr.strip()}"
        ) from exc
    finally:
        # _safe_remove(aln_path)
        _safe_remove(dnd_path)
        _safe_remove(dnd_path_cwd)


# ---------------------------------------------------------------------------
# Overlap helper
# ---------------------------------------------------------------------------
def _ranges_overlap(s1: int, e1: int, s2: int, e2: int) -> bool:
    """Return True if closed intervals [s1, e1] and [s2, e2] overlap."""
    return s1 <= e2 and s2 <= e1


# ---------------------------------------------------------------------------
# BlastParser — reciprocal best-hit evaluation
# ---------------------------------------------------------------------------
class BlastParser:
    """
    Evaluate whether the best reverse-BLAST hit overlaps the original
    reference miRNA, implementing the reciprocal best-hit criterion.

    Parameters
    ----------
    mirna : RNA object
        Reference miRNA with ``.start``, ``.end``, ``.chromosome``,
        ``.strand``, and ``.seq`` attributes.
    blasthits : list
        Parsed BLAST output rows (list of lists).
    msl : float
        Minimum sequence length ratio (currently unused here but kept
        for API compatibility with the caller).
    """

    def __init__(self, mirna, blasthits: List[List[str]], msl: float):
        self.start: int = mirna.start
        self.end: int = mirna.end
        self.chromosome: str = mirna.chromosome
        self.strand: str = mirna.strand
        self.refseq: str = mirna.seq
        self.blasthits = blasthits
        self.msl = msl

    def evaluate_besthit(self) -> bool:
        """
        Check whether any top-scoring re-BLAST hit overlaps the reference
        miRNA on the same chromosome.

        Only hits tied for the best bitscore are considered.

        Returns
        -------
        bool
            True if at least one top hit overlaps the reference miRNA.
        """
        if not self.blasthits:
            logger.info("Rejecting: No reciprocal BLAST hit found")
            return False

        top_score = float(self.blasthits[0][4])

        for hit in self.blasthits:
            if not hit:
                continue
            # Stop once we move past the top-scoring tier
            if float(hit[4]) < top_score:
                break

            sseqid = hit[0]
            sstrand = hit[3]

            if sstrand == "plus":
                sstart = int(hit[1])
                send = int(hit[2])
            elif sstrand == "minus":
                sstart = int(hit[2])
                send = int(hit[1])
            else:
                raise ValueError(
                    f"re-BLAST hit on unexpected strand: '{sstrand}'"
                )

            # Must be on the same contig
            if sseqid != self.chromosome:
                logger.info(
                    "Rejecting hit: contig mismatch — expected %s, found %s",
                    self.chromosome, sseqid,
                )
                continue

            logger.info(
                "miRNA: [%d, %d]  BLAST hit: [%d, %d]",
                self.start, self.end, sstart, send,
            )

            if _ranges_overlap(self.start, self.end, sstart, send):
                return True

            logger.info("Rejecting hit: no overlap with reference miRNA")

        return False

    def check_coortholog_ref(self, candidate_seq: str, out: str) -> bool:
        """
        When the best re-BLAST hit is not the reference miRNA itself,
        determine whether it is a co-ortholog by comparing pairwise
        distances: if the best BLAST hit is closer to the reference
        than to the candidate, accept the candidate.

        Parameters
        ----------
        candidate_seq : str
            Nucleotide sequence of the candidate ortholog.
        out : str
            Directory for temporary files.

        Returns
        -------
        bool
        """
        pid = os.getpid()
        tmp_out = os.path.join(out, f"{pid}.fa")

        best_hit_seq = self.blasthits[0][12].replace("-", "")

        with open(tmp_out, "w") as fh:
            fh.write(f">candidate\n{candidate_seq}\n")
            fh.write(f">reference\n{self.refseq}\n")
            fh.write(f">best_hit\n{best_hit_seq}\n")

        dm = calculate_distance_matrix(tmp_out)

        dist_hit_candidate = dm["best_hit", "candidate"]
        dist_ref_hit = dm["best_hit", "reference"]

        if dist_ref_hit < dist_hit_candidate:
            logger.info(
                "Co-ortholog check: dist(query, hit)=%.4f, "
                "dist(hit, ref)=%.4f → Accepting",
                dist_hit_candidate, dist_ref_hit,
            )
            return True

        logger.info(
            "Co-ortholog check: dist(query, hit)=%.4f, "
            "dist(hit, ref)=%.4f → Rejecting",
            dist_hit_candidate, dist_ref_hit,
        )
        return False


# ---------------------------------------------------------------------------
# ReBlastParser — co-ortholog verification among multiple candidates
# ---------------------------------------------------------------------------
class ReBlastParser:
    """
    Given multiple re-BLAST-verified candidates, determine which are
    true co-orthologs by comparing their sequence distances to the
    reference.

    Parameters
    ----------
    mirna : RNA object
        Reference miRNA (must have ``.seq``).
    reblast_dict : dict
        ``{candidate_id: sequence}`` for all accepted re-BLAST hits.
    """

    def __init__(self, mirna, reblast_dict: Dict[str, str]):
        self.refseq: str = mirna.seq
        self.hits = reblast_dict
        self.best_candidate: str = next(iter(reblast_dict))

    def verify_coorthologs(self, out: str) -> Dict[str, str]:
        """
        Align all candidates with the reference and accept those whose
        distance to the best candidate is smaller than the best
        candidate's distance to the reference.

        The best candidate is always included in the output.

        Parameters
        ----------
        out : str
            Directory for temporary files.

        Returns
        -------
        dict
            ``{candidate_id: sequence}`` for accepted (co-)orthologs.
        """
        pid = os.getpid()
        tmp_out = os.path.join(out, f"{pid}.fa")

        with open(tmp_out, "w") as fh:
            fh.write(f">reference\n{self.refseq}\n")
            for cand_id, seq in self.hits.items():
                fh.write(f">{cand_id}\n{seq}\n")

        dm = calculate_distance_matrix(tmp_out)

        out_dict: Dict[str, str] = {}
        dist_best_ref = dm[self.best_candidate, "reference"]

        for candidate, seq in self.hits.items():
            if candidate == self.best_candidate:
                out_dict[candidate] = seq
                continue

            dist_candidates = dm[self.best_candidate, candidate]
            if dist_candidates < dist_best_ref:
                logger.info(
                    "Co-ortholog detected: dist(%s, %s)=%.4f < "
                    "dist(%s, reference)=%.4f",
                    self.best_candidate, candidate, dist_candidates,
                    self.best_candidate, dist_best_ref,
                )
                out_dict[candidate] = seq

        return out_dict