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

# Extract genomic sequences for CM/pHMM search hits.

import logging
from typing import Dict, Iterable, Tuple

import pyfaidx


logger = logging.getLogger("ncortho")


class GenomeParser:
    """
    Extract genomic sequences for a list of CM/pHMM search hits.

    Parameters
    ----------
    genpath : str
        Path to the genome FASTA file.
    hitlist : iterable
        Each element is a tuple/list with at least 5 fields:
        ``(candidate_id, chrom, start, end, strand, ...)``.
        *start* and *end* are 1-based inclusive coordinates.
        For minus-strand hits, *start* < *end* is expected (the caller
        must normalise this beforehand).
    """

    def __init__(self, genpath: str, hitlist: Iterable):
        self.genpath = genpath
        self.hitlist = list(hitlist)
        self.genome = pyfaidx.Fasta(genpath)

    def extract_sequences(self) -> Dict[str, str]:
        """
        Return ``{candidate_id: sequence}`` for every hit.

        Plus-strand hits are returned as-is; minus-strand hits are
        reverse-complemented.

        Raises
        ------
        ValueError
            If a hit has an unrecognised strand value.
        """
        seq_dict: Dict[str, str] = {}

        for hit in self.hitlist:
            cand_id = hit[0]
            chrom = hit[1]
            start = int(hit[2])
            end = int(hit[3])
            strand = hit[4]

            # Validate coordinates
            if start > end:
                logger.warning(
                    "Hit %s has start (%d) > end (%d); swapping.",
                    cand_id, start, end,
                )
                start, end = end, start

            # pyfaidx uses 0-based half-open slicing, so subtract 1 from
            # the 1-based inclusive start.
            region = self.genome[chrom][start - 1 : end]

            if strand == "+":
                seq = region.seq
            elif strand == "-":
                seq = region.reverse.complement.seq
            else:
                raise ValueError(
                    f"Unexpected strand '{strand}' for hit {cand_id}"
                )

            seq_dict[cand_id] = seq

        return seq_dict