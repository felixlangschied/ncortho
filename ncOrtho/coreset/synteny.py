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

# Identify syntenic regions in core species and extract their sequences
# for downstream alignment and model construction.

import logging
import os
from typing import Any, Dict, List, Optional, Tuple

import pyfaidx

try:
    from coreset_utils import parse_annotation
except ModuleNotFoundError:
    from ncOrtho.coreset.coreset_utils import parse_annotation


logger = logging.getLogger("ncortho")


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------
def _vprint(msg: str, verbose: bool) -> None:
    if verbose:
        print(msg, flush=True)


def _read_genome(path: str, coretax: str, outdir: str) -> pyfaidx.Fasta:
    """
    Create a symlink to the core genome (for pyfaidx index write
    permissions) and return a ``pyfaidx.Fasta`` handle.
    """
    core_gen_dir = os.path.join(outdir, "core_genomes")
    os.makedirs(core_gen_dir, exist_ok=True)
    slink = os.path.join(core_gen_dir, coretax)
    try:
        os.symlink(path, slink)
    except FileExistsError:
        pass
    return pyfaidx.Fasta(slink)


# ---------------------------------------------------------------------------
# Synteny verification
# ---------------------------------------------------------------------------
def _synteny_check(
    left: List[str],
    right: List[str],
    orthodict: Dict[str, Any],
    mgi: int,
    verbose: bool,
) -> List[Tuple[str, int, int]]:
    """
    For each pair of (left_ortho, right_ortho), check whether they are
    on the same chromosome and within *mgi* gene positions of each
    other.  If so, return the intervening genomic region.

    Parameters
    ----------
    left : list[str]
        Ortholog gene IDs for the left anchor.
    right : list[str]
        Ortholog gene IDs for the right anchor.
    orthodict : dict
        Core-species annotation dict (``gene_id → (chrom, pos)`` and
        ``chrom → {pos: (gene, start, end, strand)}``).
    mgi : int
        Maximum gene insertions (distance in gene positions).
    verbose : bool
        Print diagnostic output.

    Returns
    -------
    list of (chrom, region_start, region_end)
    """
    regions: List[Tuple[str, int, int]] = []

    for left_ortho in left:
        if left_ortho not in orthodict:
            continue
        left_chrom, left_pos = orthodict[left_ortho]

        for right_ortho in right:
            if right_ortho not in orthodict:
                continue
            right_chrom, right_pos = orthodict[right_ortho]

            distance = abs(left_pos - right_pos)
            if left_chrom != right_chrom or distance > mgi:
                continue

            left_info = orthodict[left_chrom][left_pos]
            right_info = orthodict[left_chrom][right_pos]
            # info tuple: (gene_id, start, end, strand)
            left_start, left_end = int(left_info[1]), int(left_info[2])
            right_start, right_end = int(right_info[1]), int(right_info[2])

            # Extract the region *between* the two anchor genes
            if left_pos < right_pos:
                region_start = left_end
                region_end = right_start
            else:
                region_start = right_end
                region_end = left_start

            if region_end <= region_start:
                _vprint(
                    f"Skipping overlapping/adjacent anchors "
                    f"{left_ortho}–{right_ortho} on {left_chrom}",
                    verbose,
                )
                continue

            regions.append((left_chrom, region_start, region_end))

    return regions


# ---------------------------------------------------------------------------
# Main entry point
# ---------------------------------------------------------------------------
def analyze_synteny(
    core_d: Dict[str, Dict[str, str]],
    mirna_pos: Dict[str, Tuple],
    out: str,
    idtype: str,
    mgi: int,
    verbose: bool,
) -> Dict[str, List[str]]:
    """
    For every miRNA × core species pair, identify syntenic regions and
    extract their genomic sequences.

    Parameters
    ----------
    core_d : dict
        ``{species: {"genome": path, "annotation": path, ...}}``.
    mirna_pos : dict
        ``{mirid: (synteny_type, {taxon: (left_orthos, right_orthos)})}``.
    out : str
        Base output / tmp directory.
    idtype : str
        Gene-ID type for annotation parsing.
    mgi : int
        Maximum gene insertions allowed between anchors.
    verbose : bool
        Print diagnostic output.

    Returns
    -------
    dict
        ``{mirid: [">taxon_0\\n", "sequence\\n", ...]}``.
    """
    collector: Dict[str, List[str]] = {}

    for taxon, spec_data in core_d.items():
        print(f"# {taxon}", flush=True)
        _vprint(f"# Parsing annotation file for {taxon}", verbose)
        core_anno = parse_annotation(spec_data["annotation"], idtype)

        _vprint("# Loading genome file", verbose)
        genome = _read_genome(spec_data["genome"], taxon, out)
        _vprint("# Done", verbose)

        for mirid, (style, ortho_per_taxon) in mirna_pos.items():
            if mirid not in collector:
                collector[mirid] = []

            if taxon not in ortho_per_taxon:
                _vprint(f"No core orthologs found for {mirid} in {taxon}", verbose)
                continue

            left_orthos, right_orthos = ortho_per_taxon[taxon]
            found = False

            # --- "inside" / "opposite": left and right anchors are the
            # same gene — extract the ortholog's own extent. ----------
            if left_orthos == right_orthos:
                for idx, ortholog in enumerate(left_orthos):
                    if ortholog not in core_anno:
                        continue
                    ortho_chrom, ortho_pos = core_anno[ortholog]

                    ortho_info = core_anno[ortho_chrom][ortho_pos]
                    ortho_start = int(ortho_info[1])
                    ortho_end = int(ortho_info[2])

                    seq = genome[ortho_chrom][ortho_start - 1 : ortho_end].seq
                    collector[mirid].append(f">{taxon}_{idx}\n")
                    collector[mirid].append(f"{seq}\n")
                    found = True

            # --- "in-between": extract region between left and right
            # anchor genes. ----------------------------------------
            else:
                syn_regions = _synteny_check(
                    left_orthos, right_orthos, core_anno, mgi, verbose,
                )
                for count, (syn_chrom, syn_start, syn_end) in enumerate(syn_regions):
                    seq = genome[syn_chrom][syn_start - 1 : syn_end].seq
                    collector[mirid].append(f">{taxon}_{count}\n")
                    collector[mirid].append(f"{seq}\n")
                    found = True

            if found:
                _vprint(f"Synteny fulfilled for {mirid} in {taxon}", verbose)
            else:
                _vprint(f"No syntenic region found for {mirid} in {taxon}", verbose)

    return collector