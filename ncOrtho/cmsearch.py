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

# CM / pHMM search pipeline:
#   heuristic BLAST -> candidate region extraction -> model search -> parsing -> deduplication

import logging
import os
import subprocess as sp
from typing import Dict, List, Optional, Tuple

import pyfaidx


logger = logging.getLogger("ncortho")

# Default flanking region (bp) added around BLAST hits when extracting
# candidate sequences for the CM/pHMM search.
DEFAULT_EXTRA_REGION = 1000


# ---------------------------------------------------------------------------
# BLAST-based heuristic candidate search
# ---------------------------------------------------------------------------
def _run_candidate_blast(
    blastdb: str,
    cpu: int,
    evalue: float,
    sequence: str,
    min_hit_length: float,
    task: str = "blastn",
) -> List[List[str]]:
    """
    BLAST *sequence* against *blastdb* and return hits whose alignment
    length is at least *min_hit_length*.

    Returns a list of hits, each being ``[sseqid, sstart, send, sstrand, length]``.
    """
    blast_cmd = [
        "blastn",
        "-task", task,
        "-db", blastdb,
        "-num_threads", str(cpu),
        "-evalue", str(evalue),
        "-outfmt", "6 sseqid sstart send sstrand length",
    ]
    result = sp.run(
        blast_cmd,
        input=sequence,
        capture_output=True,
        text=True,
        check=False,
    )
    if result.returncode != 0:
        raise sp.SubprocessError(
            f"BLASTn failed (exit {result.returncode}): {result.stderr.strip()}"
        )

    hit_list = []
    for line in result.stdout.strip().split("\n"):
        if not line:
            continue
        fields = line.split()
        hit_length = float(fields[-1])
        if hit_length >= min_hit_length:
            hit_list.append(fields)
    return hit_list


def _extract_candidate_regions(
    candidate_list: List[List[str]],
    query: str,
    extra_region: int = DEFAULT_EXTRA_REGION,
) -> List[str]:
    """
    For each BLAST hit, extract the surrounding genomic region (hit ±
    *extra_region* bp) and return a list of FASTA-formatted strings.

    The FASTA header encodes the original coordinates so that CM hits
    can later be mapped back to the genome.
    """
    genome = pyfaidx.Fasta(query)
    regions: List[str] = []

    for hit in candidate_list:
        chrom, start_s, end_s, strand, _length = hit
        strand = strand.replace("plus", "+").replace("minus", "-")

        start = int(start_s)
        end = int(end_s)
        if strand == "-":
            start, end = end, start

        hit_at_start = False
        hit_at_end = False

        chrom_length = genome[chrom][-1].end

        if start > extra_region:
            n_start = start - extra_region
        else:
            n_start = 0
            hit_at_start = True

        n_end = end + extra_region
        if n_end > chrom_length:
            n_end = chrom_length
            hit_at_end = True

        if strand == "+":
            sequence = genome[chrom][n_start:n_end].seq
        elif strand == "-":
            sequence = genome[chrom][n_start:n_end].reverse.complement.seq
        else:
            raise ValueError(f"Unknown strand: {strand}")

        header = f">{chrom}|{start}|{end}|{strand}|{hit_at_start}|{hit_at_end}"
        regions.append(f"{header}\n{sequence}\n")

    return regions


def heuristic_search(
    blastdb: str,
    rna,
    query: str,
    heuristic: Tuple,
    cpu: int,
    out: str,
) -> str:
    """
    Run a preliminary BLAST search of the reference miRNA against the
    query genome to find candidate regions for the downstream CM search.

    Returns the path to a FASTA file with candidate regions, or an empty
    string if no candidates pass the length cutoff.
    """
    blast_len_cut = len(rna.seq) * heuristic[2]

    candidate_list = _run_candidate_blast(
        blastdb, cpu, heuristic[1], rna.seq, blast_len_cut,
    )
    if not candidate_list:
        return ""

    candidate_regions = _extract_candidate_regions(candidate_list, query)
    fasta = os.path.join(out, f"candidate_regions_{rna.name}.fa")
    with open(fasta, "w") as fh:
        for line in candidate_regions:
            fh.write(line)

    return fasta


# ---------------------------------------------------------------------------
# CM / pHMM model search
# ---------------------------------------------------------------------------
def _parse_phmm_output(tblout: str) -> List[List[str]]:
    """
    Parse nhmmer ``--tblout`` output into the same column layout that
    ``_parse_cm_output`` returns so downstream code can be agnostic.
    """
    results = []
    with open(tblout) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            data = line.strip().split()
            info = data[0]
            start, end = data[6:8]
            strand, evalue, score = data[11:14]
            # Pad with placeholder columns to match CM output layout:
            # [rank, type, evalue, score, bias, info, start, end, strand,
            #  trunc, gc, noncannuc, ...]
            results.append(
                ["_", "_", evalue, score, "_", info, start, end, strand, "_", "_", "_", "_"]
            )
    return results


def _parse_cm_output(stdout: str) -> List[List[str]]:
    """Parse cmsearch stdout (hit table) into a list of hit rows."""
    results = []
    for line in stdout.split("\n"):
        if line.startswith("  ("):
            if "No hits detected that satisfy reporting thresholds" in line:
                return []
            results.append(line.strip().split())
    return results


def perform_model_search(
    models: str,
    rna,
    target_fasta: str,
    output: str,
    cm_cutoff: float,
    cpu: int,
    phmm: bool,
) -> List[List[str]]:
    """
    Run cmsearch (or nhmmer for pHMM) and return the raw parsed hit list.
    """
    cutoff = rna.bit * cm_cutoff

    if phmm:
        model_file = os.path.join(models, f"{rna.name}.phmm")
        tblout = os.path.join(output, f"{rna.name}_phmm.out")
        cmd = [
            "nhmmer",
            "--cpu", str(cpu),
            "--tblout", tblout,
            "-T", str(cutoff),
            "--incT", str(cutoff),
            model_file,
            target_fasta,
        ]
        result = sp.run(cmd, capture_output=True, text=True, check=False)
        if result.returncode != 0:
            raise sp.SubprocessError(
                f"nhmmer failed (exit {result.returncode}): {result.stderr.strip()}"
            )
        results = _parse_phmm_output(tblout)
        os.remove(tblout)
    else:
        model_file = os.path.join(models, f"{rna.name}.cm")
        cmd = [
            "cmsearch",
            "--cpu", str(cpu),
            "--noali",
            "-T", str(cutoff),
            "--incT", str(cutoff),
            model_file,
            target_fasta,
        ]
        result = sp.run(cmd, capture_output=True, text=True, check=False)
        if result.returncode != 0:
            raise sp.SubprocessError(
                f"cmsearch failed (exit {result.returncode}): {result.stderr.strip()}"
            )
        results = _parse_cm_output(result.stdout)

    return results


# ---------------------------------------------------------------------------
# CM result coordinate parsing
# ---------------------------------------------------------------------------
def parse_cmsearch(cm_results: List[List[str]]) -> List[List]:
    """
    Parse CM results from a direct (non-heuristic) search into
    ``[chrom, start, end, strand, score]`` rows.
    """
    parsed = []
    for data in cm_results:
        score = float(data[3])
        chrom = data[5]
        start, end = int(data[6]), int(data[7])
        strand = data[8]
        parsed.append([chrom, start, end, strand, score])
    return parsed


def parse_cmsearch_for_heuristic(
    cm_results: List[List[str]],
    extra_region: int = DEFAULT_EXTRA_REGION,
) -> List[List]:
    """
    Map CM hits from candidate-region-local coordinates back to
    genome-global coordinates.

    The FASTA header of each candidate region encodes the original BLAST
    hit position, boundary flags, and strand as:
    ``chrom|blast_start|blast_end|strand|hit_at_start|hit_at_end``
    """
    parsed = []

    for data in cm_results:
        # Decode the encoded FASTA header
        header_parts = data[5].split("|")
        blast_chrom = header_parts[0]
        blast_start = int(header_parts[1])
        blast_end = int(header_parts[2])
        blast_strand = header_parts[3]
        hit_at_start = header_parts[4] == "True"
        hit_at_end = header_parts[5] == "True"

        cm_start = int(data[6])
        cm_end = int(data[7])
        cm_strand = data[8]

        # Reverse CM hits should not occur since BLAST hits on the minus
        # strand were already reverse-complemented before the CM search.
        if cm_strand == "-":
            continue

        # Convert candidate-region-local coordinates to genome-global.
        if hit_at_start:
            # Region started at position 0 (no left flank), so CM coords
            # are already genome-global relative to the BLAST hit start.
            hit_start = cm_start
            hit_end = cm_end
        else:
            # Normal case: the region starts at (blast_start - extra_region),
            # so offset the CM positions accordingly.
            region_origin = blast_start - extra_region
            hit_start = region_origin + cm_start - 1
            hit_end = region_origin + cm_end - 1

        # NOTE: When hit_at_end is True, the right flank was truncated,
        # but the origin calculation is the same as the normal case —
        # only the left boundary matters for the offset.

        parsed.append([blast_chrom, hit_start, hit_end, blast_strand, float(data[3])])

    return parsed


# ---------------------------------------------------------------------------
# Deduplication and cutoff filtering
# ---------------------------------------------------------------------------
def _ranges_overlap(s1: int, e1: int, s2: int, e2: int) -> bool:
    """Return True if [s1, e1] and [s2, e2] overlap (inclusive)."""
    return s1 <= e2 and s2 <= e1


def remove_duplicates_apply_cutoffs(
    cms: List[List],
    cm_cutoff: float,
    length_cutoff: float,
    rna,
) -> Dict[str, Tuple]:
    """
    Filter CM hits by score and length, then remove overlapping
    duplicates on opposite strands (keeping the higher-scoring one).

    Parameters
    ----------
    cms : list
        CM hits as ``[chrom, start, end, strand, score]``.
    cm_cutoff : float
        Minimum bitscore ratio relative to the reference score.
    length_cutoff : float
        Minimum hit length as a fraction of the reference miRNA length.
    rna : RNA object
        Must have ``.bit`` (reference bitscore) and ``.seq`` attributes.

    Returns
    -------
    dict
        ``{candidate_id: (candidate_id, chrom, start, end, strand, score)}``
    """
    min_score = rna.bit * cm_cutoff
    min_length = len(rna.seq) * length_cutoff

    hits_dict: Dict[str, Tuple] = {}
    # Group hits by chromosome for overlap detection
    chromo_dict: Dict[str, List[Tuple]] = {}

    for candidate_nr, hit in enumerate(cms, 1):
        h_chrom, h_start, h_end, h_strand, h_score = hit

        # Ensure start <= end (blastparser expects this)
        if h_start > h_end:
            h_start, h_end = h_end, h_start
        length = h_end - h_start

        # BUG-FIX: Compare against the *reference* miRNA length, not
        # against the hit's own length (which was always True before).
        if length < min_length:
            continue
        if h_score < min_score:
            continue

        candidate = f"{rna.name}_c{candidate_nr}"
        entry = (candidate, h_chrom, h_start, h_end, h_strand, h_score)
        hits_dict[candidate] = entry

        # BUG-FIX: Append to the per-chromosome list instead of
        # overwriting it. The original code always reset the list to a
        # single element, so duplicate detection across multiple hits on
        # the same chromosome never worked.
        if h_chrom not in chromo_dict:
            chromo_dict[h_chrom] = []
        chromo_dict[h_chrom].append(entry)

    # Remove overlapping hits on opposite strands, keeping the higher score.
    to_remove = set()
    for chrom, chrom_hits in chromo_dict.items():
        n = len(chrom_hits)
        for i in range(n):
            cand_i = chrom_hits[i]
            _, _, start_i, end_i, strand_i, score_i = cand_i
            for j in range(i + 1, n):
                cand_j = chrom_hits[j]
                _, _, start_j, end_j, strand_j, score_j = cand_j
                # Only check opposite-strand overlaps
                if strand_i == strand_j:
                    continue
                if _ranges_overlap(start_i, end_i, start_j, end_j):
                    if score_i >= score_j:
                        to_remove.add(cand_j[0])
                    else:
                        to_remove.add(cand_i[0])

    for cand_id in to_remove:
        hits_dict.pop(cand_id, None)

    return hits_dict


# ---------------------------------------------------------------------------
# Top-level search dispatcher
# ---------------------------------------------------------------------------
def model_search(
    rna,
    cm_cutoff: float,
    cpu: int,
    msl: float,
    models: str,
    query: str,
    blastdb: str,
    out: str,
    cleanup: bool,
    heuristic_col: Tuple,
    phmm: bool,
) -> Tuple[Optional[Dict], str]:
    """
    Run the full model-search pipeline for a single miRNA.

    Steps:
      1. (Optional) Heuristic BLAST to narrow candidate regions.
      2. CM or pHMM search in candidate regions (or full genome).
      3. Parse and remap coordinates.
      4. Deduplicate and apply cutoffs.

    Returns
    -------
    (hits_dict, status_message) where *hits_dict* is ``None`` on failure.
    """
    # Check that the model file exists
    if phmm:
        model_file = os.path.join(models, f"{rna.name}.phmm")
    else:
        model_file = os.path.join(models, f"{rna.name}.cm")

    if not os.path.isfile(model_file):
        return None, "No model found"

    heuristic, _heur_eval, _heur_len, sensitive_heuristic = heuristic_col
    candidate_region_fasta = ""

    if heuristic:
        candidate_region_fasta = heuristic_search(
            blastdb, rna, query, heuristic_col, cpu, out,
        )

        if not candidate_region_fasta and sensitive_heuristic:
            # No candidates above cutoff; fall back to full-genome search.
            cm_results = perform_model_search(
                models, rna, query, out, cm_cutoff, cpu, phmm,
            )
            cm_results = parse_cmsearch(cm_results)
        elif not candidate_region_fasta:
            return None, "No candidate region found with BLASTn"
        else:
            cm_results = perform_model_search(
                models, rna, candidate_region_fasta, out, cm_cutoff, cpu, phmm,
            )
            cm_results = parse_cmsearch_for_heuristic(cm_results)

        # BUG-FIX: Only remove the file if it actually exists (it is an
        # empty string when the sensitive_heuristic fallback was used).
        if cleanup and candidate_region_fasta and os.path.isfile(candidate_region_fasta):
            os.remove(candidate_region_fasta)
    else:
        cm_results = perform_model_search(
            models, rna, query, out, cm_cutoff, cpu, phmm,
        )
        cm_results = parse_cmsearch(cm_results)

    if not cm_results:
        return None, "No CMsearch results in candidate regions or query genome"

    parsed = remove_duplicates_apply_cutoffs(cm_results, cm_cutoff, msl, rna)
    if not parsed:
        return None, "No CMsearch results above threshold"

    return parsed, "Success"