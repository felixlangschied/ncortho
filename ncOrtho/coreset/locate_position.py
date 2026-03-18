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

# Determine the genomic context (position category) of each reference
# miRNA relative to protein-coding genes, and identify the corresponding
# syntenic anchor orthologs in core species.

import logging
from typing import Any, Dict, List, Optional, Tuple

logger = logging.getLogger("ncortho")


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------
def _vprint(msg: str, verbose: bool) -> None:
    if verbose:
        print(msg, flush=True)


# ---------------------------------------------------------------------------
# Boundary check
# ---------------------------------------------------------------------------
def _no_synteny_possible(
    mirna: str,
    chrom: str,
    start: int,
    end: int,
    reference: Dict[str, Any],
) -> bool:
    """
    Return True if synteny analysis is impossible for this miRNA because
    it falls outside the annotated gene range on its contig.

    Covers three cases:
      1. No protein-coding genes on the contig at all.
      2. miRNA is left of the first annotated gene.
      3. miRNA is right of the last annotated gene.
    """
    if chrom not in reference:
        logger.warning(
            'No protein-coding genes found on contig "%s". '
            "Verify that contig IDs in the miRNA input match the reference annotation.",
            chrom,
        )
        return True

    chrom_genes = reference[chrom]
    first_gene_info = chrom_genes[1]
    first_gene_start = int(first_gene_info[1])
    last_gene_info = chrom_genes[len(chrom_genes)]
    last_gene_end = int(last_gene_info[2])

    if end < first_gene_start:
        logger.info(
            "No left neighbor of %s (start of contig %s). "
            "Right neighbor: %s",
            mirna, chrom, first_gene_info[0],
        )
        return True

    if start > last_gene_end:
        logger.info(
            "No right neighbor of %s (end of contig %s). "
            "Left neighbor: %s",
            mirna, chrom, last_gene_info[0],
        )
        return True

    return False


# ---------------------------------------------------------------------------
# Ortholog neighbour search
# ---------------------------------------------------------------------------
def _find_ortho(
    genename: str,
    position: int,
    chrom_genes: Dict[int, Tuple],
    orthologs: Dict[str, List[str]],
    max_distance: int,
    direction: str,
    verbose: bool,
) -> Optional[List[str]]:
    """
    Find the ortholog(s) for *genename* in a core species, or — if no
    ortholog exists — try up to *max_distance* neighbouring genes in
    the given *direction* (``"left"`` or ``"right"``).
    """
    if genename in orthologs:
        _vprint(
            f'{" ".join(orthologs[genename])} identified as ortholog(s) to {genename}',
            verbose,
        )
        return orthologs[genename]

    num_genes = len(chrom_genes)
    for offset in range(1, max_distance):
        if direction == "left":
            idx = position - offset
        else:
            idx = position + offset

        if idx < 1 or idx > num_genes:
            continue

        next_gene = chrom_genes[idx][0]
        if next_gene in orthologs:
            _vprint(
                f'No ortholog(s) to {genename}; '
                f'{" ".join(orthologs[next_gene])} identified as ortholog(s) '
                f"of {next_gene} ({offset} gene(s) away)",
                verbose,
            )
            return orthologs[next_gene]

    return None


def _neighbor_search(
    leftgene: str,
    rightgene: str,
    core_dict: Dict[str, Dict[str, List[str]]],
    gene_position: int,
    chrom_genes: Dict[int, Tuple],
    max_distance: int,
    verbose: bool,
) -> Dict[str, Tuple[List[str], List[str]]]:
    """
    For each core species, find orthologs of the left and right anchor
    genes (or their nearest annotated neighbours within *max_distance*).

    A species is only included in the result if *both* anchors can be
    resolved.
    """
    result: Dict[str, Tuple[List[str], List[str]]] = {}

    for core_taxon, orthodict in core_dict.items():
        left_orthos = _find_ortho(
            leftgene, gene_position, chrom_genes,
            orthodict, max_distance, "left", verbose,
        )
        right_orthos = _find_ortho(
            rightgene, gene_position + 1, chrom_genes,
            orthodict, max_distance, "right", verbose,
        )
        if left_orthos and right_orthos:
            result[core_taxon] = (left_orthos, right_orthos)

    return result


# ---------------------------------------------------------------------------
# Main categorisation
# ---------------------------------------------------------------------------
def categorize_mirna_position(
    mirna: str,
    mirna_chrom: str,
    mirna_start: int,
    mirna_end: int,
    mirna_strand: str,
    reference: Dict[str, Any],
    all_orthologs: Dict[str, Dict[str, List[str]]],
    no_add_orthos: int,
    verbose: bool,
) -> Tuple[Optional[str], Optional[Dict]]:
    """
    Determine the genomic context of a reference miRNA and identify
    syntenic anchor orthologs in core species.

    Position categories:
      - ``"inside"``:     miRNA is within a gene on the same strand.
      - ``"opposite"``:   miRNA overlaps a gene on the opposite strand.
      - ``"in-between"``: miRNA lies between two adjacent genes.

    Parameters
    ----------
    reference : dict
        Annotation dict as returned by the parsers in ``coreset_utils``.
    all_orthologs : dict
        ``{core_taxon: {ref_gene: [ortho, ...], ...}, ...}``.
    no_add_orthos : int
        Number of additional genes to consider as alternative anchors.
    verbose : bool
        Print diagnostic messages.

    Returns
    -------
    (synteny_type, core_orthologs)
        ``(None, None)`` if the position cannot be resolved.
    """
    if _no_synteny_possible(mirna, mirna_chrom, mirna_start, mirna_end, reference):
        return None, None

    chrom_genes = reference[mirna_chrom]

    for position, geneinfo in chrom_genes.items():
        if not isinstance(position, int):
            continue

        gene, gene_start, gene_end, gene_strand = geneinfo

        # Case: miRNA is entirely within a gene
        if mirna_start >= gene_start and mirna_end <= gene_end:
            if mirna_strand == gene_strand:
                syntype = "inside"
                _vprint(f"{mirna} is located inside gene {gene}", verbose)
            else:
                syntype = "opposite"
                _vprint(f"{mirna} is located opposite of gene {gene}", verbose)

            core_orthos = _neighbor_search(
                gene, gene, all_orthologs, position,
                chrom_genes, no_add_orthos, verbose,
            )
            return syntype, core_orthos

        # Case: miRNA is between this gene and the next
        if position + 1 not in chrom_genes:
            continue

        right_info = chrom_genes[position + 1]
        rightneighbor, _rn_start, rn_end, _rn_strand = right_info

        if gene_start <= mirna_start and rn_end >= mirna_end:
            _vprint(
                f"{gene} is the left neighbor of {mirna}, "
                f"{rightneighbor} is the right neighbor",
                verbose,
            )
            core_orthos = _neighbor_search(
                gene, rightneighbor, all_orthologs, position,
                chrom_genes, no_add_orthos, verbose,
            )
            return "in-between", core_orthos

    _vprint(f"Unable to resolve synteny for {mirna}", verbose)
    return None, None