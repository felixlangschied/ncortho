'''
ncOrtho submodule
Extract (pre-miRNA) sequence from query genome according to cmsearch
coordinates
Makes use of pyfaidx module for efficient sequence extraction from
FASTA files
https://pypi.org/project/pyfaidx/
TODO: include license, author details etc
'''

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
                    self.gene_dict[hit[1]][int(hit[3])-1:int(hit[2])]
                    .reverse.complement.seq
                )
            seq_dict[hit[0]] = seq
        return seq_dict
