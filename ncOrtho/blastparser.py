import subprocess as sp
from Bio import AlignIO
from Bio.Phylo.TreeConstruction import DistanceCalculator
import os

'''
ncOrtho submodule
TODO: include license, author details etc

### TODO: calculate length of overlap to see if it is relevant/significant
### also the bit score should be compared
### add a check for the presence of the mature miRNA
'''

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
    def __init__(self, mirna, blastpath, msl):
        self.start = mirna.start
        self.end = mirna.end
        self.chromosome = mirna.chromosome
        self.strand = mirna.strand
        self.refseq = mirna.pre
        del mirna
        #self.blastpath = blastpath
        with open(blastpath, 'r') as blastfile:
            self.blasthits = [line.strip().split() for line in blastfile]
        self.msl = msl
        #self.blasthits = []
        #self.top_score = 100
        #self.top_score = blasthits[0][6]

    def evaluate_besthit(self,):
        # with open(self.blastpath) as blastfile:
        #     blasthits = [line.strip().split() for line in blastfile]
        # If no BLAST hit was found, the search failed by default.
        if not self.blasthits:
            print('No reciprocal BLAST hit found')
            return False
        # Otherwise, check if the best hit and the reference miRNA overlap.
        else:
            # Gather coordinates of best BLAST hit.
            tophit = self.blasthits[0]
            sseqid = tophit[1]
            sstart = int(tophit[8])
            send = int(tophit[9])
            # del blasthits
            # Sequences must be on the same contig, otherwise overlap can be
            # ruled out instantaneously
            if not sseqid == self.chromosome:
                return False
            # Contigs match, so overlap is possible.
            else:
                # first within second
                if (
                    (sstart <= self.start and self.start <= send)
                    or (sstart <= self.end and self.end <= send)
                ):
                    return True
                # second within first
                elif (
                    (self.start <= sstart and sstart <= self.end)
                    or (self.start <= send and send <= self.end)
                ):
                    return True
                # No overlap
                else:
                    return False

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
