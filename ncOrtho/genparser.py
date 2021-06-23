"""
ncOrtho - Targeted ortholog search for miRNAs
Copyright (C) 2021 Felix Langschied

ncOrtho is a free software: you can redistribute it and/or modify
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

import pyfaidx

class GenomeParser(object):
    
    # genpath: path to the genome file to extract from
    # hitlist: list of hits to extract the sequence for
    def __init__(self, genpath, hitlist):
        self.genpath = genpath
        self.hitlist = hitlist
        self.gene_dict = self.parse_genome()
    
    # read in the genome
    def parse_genome(self,):
        genome = pyfaidx.Fasta(self.genpath)
        return genome
    
    # extract the sequences for each significant hit from cmsearch
    def extract_sequences(self,):
        # dictionary to store the sequences
        seq_dict = {}
        for hit in self.hitlist:
            # hit on sense strand, sequence can be extracted directly
            if hit[4] == '+':
                seq = self.gene_dict[hit[1]][int(hit[2])-1:int(hit[3])].seq
            # hit on antisense strand, sequence has to be extracted
            # as reverse complement
            # note that cmsearch switches start and end positions in case
            # of antisense strand hits
            elif hit[4] == '-':
                seq = (
                    self.gene_dict[hit[1]][int(hit[2])-1:int(hit[3])]
                    .reverse.complement.seq
                )
            seq_dict[hit[0]] = seq
        return seq_dict
