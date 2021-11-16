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

import subprocess as sp
from Bio import AlignIO
from Bio.Phylo.TreeConstruction import DistanceCalculator
import os
import sys


def calculate_distance_matrix(aln_path):
    # align reBLAST hits
    aln_cmd = (
        't_coffee {0} -no_warning -quiet -type=dna '
        ' -output=fasta_aln -outfile={0}'.format(aln_path)
    )
    sp.run(aln_cmd, shell=True)
    aln = AlignIO.read(open(aln_path), 'fasta')
    calculator = DistanceCalculator('identity')
    dm = calculator.get_distance(aln)
    # delete temporary files
    os.remove(aln_path)
    curr_path = os.getcwd()
    os.remove('{}/{}.dnd'.format(curr_path, aln_path.split('/')[-1].split('.fa')[0]))
    return dm


# BlastParser object performs reverse BLAST search and reports if a candidate
# fulfills the reverse best hit criterion
class BlastParser(object):
    # init parameters:
    # mirna: Mirna object that holds the location for the reference miRNA
    # blastpath: path to the output file of the reverse BLAST search
    # msl: ncOrtho minimum sequence length threshold
    def __init__(self, mirna, blasthits, msl):
        self.start = mirna.start
        self.end = mirna.end
        self.chromosome = mirna.chromosome
        self.strand = mirna.strand
        self.refseq = mirna.pre
        del mirna
        self.blasthits = blasthits
        self.msl = msl

    def evaluate_besthit(self,):
        # with open(self.blastpath) as blastfile:
        #     blasthits = [line.strip().split() for line in blastfile]
        # If no BLAST hit was found, the search failed by default.
        outbool = False
        if not self.blasthits:
            print('Rejecting: No reciprocal BLAST hit found')
            return False
        # Otherwise, check if the best hit and the reference miRNA overlap.
        else:
            topscore = float(self.blasthits[0][4])
            for tophit in self.blasthits:
                if not tophit or float(tophit[4]) < topscore:
                    return False
                sseqid = tophit[0]
                sstrand = tophit[3]
                if sstrand == 'plus':
                    sstart = int(tophit[1])
                    send = int(tophit[2])
                elif sstrand == 'minus':
                    sstart = int(tophit[2])
                    send = int(tophit[1])
                else:
                    raise ValueError('re-BLAST on neither plus nor minus strand')

                # Sequences must be on the same contig, otherwise overlap can be
                # ruled out instantaneously
                if not sseqid == self.chromosome:
                    print(f'Rejecting: Contig/Chromosome does not match. Expected {self.chromosome} but found {sseqid}')
                # Contigs match, so overlap is possible.
                else:
                    print(f'Start miRNA: {self.start} End miRNA: {self.end}')
                    print(f'Start BLAST hit: {sstart} End BLAST hit: {send}')
                    # first within second
                    if (
                        (sstart <= self.start <= send)
                        or (sstart <= self.end <= send)
                    ):
                        return True
                    # second within first
                    elif (
                        (self.start <= sstart <= self.end)
                        or (self.start <= send <= self.end)
                    ):
                        return True
                    # No overlap
                    else:
                        print('Rejecting: No overlap between miRNA and best BLAST hit')
        return outbool

    def check_coortholog_ref(self, candidate_seq, out):
        # 1) get query sequence
        # 2) get reference sequence
        # 3) get sequence of best blast hit
        # 4) align sequences
        # 5) calculate difference
        # 6) if distance between best blast hit and reference smaller than
        # difference between best hit and candidate, accept candidate -> reBLAST hit co-ortholog instead of reference

        # get process id to ensure that the correct library is used when running ncOrtho in parallel
        pid = os.getpid()
        tmp_out = '{0}/{1}.fa'.format(out, pid)
        # extract sequences of reBLAST hits
        with open(tmp_out, 'w') as of:
            of.write('>candidate\n')
            of.write('{}\n'.format(candidate_seq))
            of.write('>reference\n')
            of.write('{}\n'.format(self.refseq))
            of.write('>best_hit\n')
            of.write(self.blasthits[0][12].replace('-', ''))

        dm = calculate_distance_matrix(tmp_out)
        distance_hit_query = dm['best_hit', 'candidate']
        distance_ref_hit = dm['best_hit', 'reference']

        if distance_ref_hit < distance_hit_query:
            print(
                "\t Distance query - BLAST hit: %6.4f, Distance blast hit - reference: %6.4f\tAccepting\n"
                %(distance_hit_query, distance_ref_hit)
            )
            return True
        else:
            print(
                "\t Distance query - BLAST hit: %6.4f, Distance blast hit - reference: %6.4f\tRejecting\n"
                %(distance_hit_query, distance_ref_hit)
            )
            return False


class ReBlastParser(object):
    def __init__(self, mirna, reblast_dict):
        self.refseq = mirna.pre
        self.hits = reblast_dict
        self.best_candidate = list(reblast_dict.keys())[0]

    def verify_coorthologs(self, out):
        out_dict = {}
        # get process id to ensure that the correct library is used when running ncOrtho in parallel
        pid = os.getpid()
        tmp_out = '{0}/{1}.fa'.format(out, pid)
        with open(tmp_out, 'w') as of:
            of.write('>reference\n')
            of.write(f'{self.refseq}\n')
            for candidate in self.hits:
                of.write(f'>{candidate}\n')
                of.write(f'{self.hits[candidate]}\n')

        dm = calculate_distance_matrix(tmp_out)
        for candidate in self.hits:
            if candidate == self.best_candidate:
                out_dict[candidate] = self.hits[candidate]
            else:
                distance_candidates = dm[self.best_candidate, candidate]
                distance_best_ref = dm[self.best_candidate, 'reference']
                if distance_candidates < distance_best_ref:
                    print(
                        f'co-ortholog detected: distance of best candidate to {candidate} {distance_candidates} '
                        f'compared to distance of best candidate to reference of {distance_best_ref}'
                    )
                    out_dict[candidate] = self.hits[candidate]
        return out_dict
