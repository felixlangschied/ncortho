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


# Modules import
import argparse
import multiprocessing as mp
import os
import subprocess as sp
import sys

# Internal ncOrtho modules
try:
    from ncOrtho.blastparser import BlastParser
    from ncOrtho.blastparser import ReBlastParser
    from ncOrtho.genparser import GenomeParser
    from ncOrtho.cmsearch import cmsearcher
    from ncOrtho.utils import check_blastdb
    from ncOrtho.utils import make_blastndb
    from ncOrtho.utils import find_refbit
except ModuleNotFoundError:
    from blastparser import BlastParser
    from blastparser import ReBlastParser
    from genparser import GenomeParser
    from cmsearch import cmsearcher
    from utils import check_blastdb
    from utils import make_blastndb
    from utils import find_refbit

###############################################################################


# Central class of microRNA objects
class Mirna(object):
    def __init__(self, name, chromosome, start, end, strand, pre, mature, bit):
        # miRNA identifier
        self.name = name
        # chromosome that the miRNA is located on
        if chromosome.startswith('chr'):
            self.chromosome = chromosome.split('chr')[1]
        else:
            self.chromosome = chromosome
        # start position of the pre-miRNA
        self.start = int(start)
        # end position of the pre-miRNA
        self.end = int(end)
        # sense (+) or anti-sense (-) strand
        self.strand = strand
        # nucleotide sequence of the pre-miRNA
        self.pre = pre.replace('U', 'T')
        # nucleotide sequence of the mature miRNA
        self.mature = mature.replace('U', 'T')
        # reference bit score that miRNA receives by its own
        # covariance model
        self.bit = bit


def mirna_maker(mirpath, cmpath, output, msl):
    """
    Parses the miRNA data input file and returns a dictionary of Mirna objects.

    Parameters
    ----------
    mirpath : STR
        Path to file with microRNA data.
    cmpath : STR
        Path to covariance models.
    output : STR
        Path for writing temporary files.
    msl : FLOAT
        Length filter.

    Returns
    -------
    mmdict : DICT
        Dictionary with Mirna objects.

    """
    
    mmdict = {} # will be the return object
    
    with open(mirpath) as mirna_file:
        mirna_data = [
            line.strip().split() for line in mirna_file
            if not line.startswith('#')
        ]

    print('# Calculating reference bit scores')
    for mirna in mirna_data:
        mirid = mirna[0]
        seq = mirna[5]
        query = '{0}/{1}.fa'.format(output, mirid)
        model = '{0}/{1}.cm'.format(cmpath, mirid)
        # Check if the covariance model even exists, otherwise skip to
        # the next miRNA.
        if not os.path.isfile(model):
            print('WARNING: No covariance model found for {}'.format(mirid))
            mmdict[mirna[0]] = ''
            continue
        # Check if the output folder exists, otherwise create it.
        if not os.path.isdir('{}'.format(output)):
            mkdir = 'mkdir -p {}'.format(output)
            sp.call(mkdir, shell=True)

        # Obtain the reference bit score for each miRNA by applying it
        # to its own covariance model.
        
        # Create a temporary FASTA file with the miRNA sequence as
        # query for external search tool cmsearch to calculate the
        # reference bit score.
        with open(query, 'w') as tmpfile:
            tmpfile.write('>{0}\n{1}'.format(mirid, seq))
        cms_command = (
            'cmsearch -E 0.01 --noali {0} {1}'
            .format(model, query)
        )
        res = sp.run(cms_command, shell=True, capture_output=True)
        if res.returncode == 0:
            top_score = find_refbit(res.stdout.decode('utf-8'))
            mirna.append(top_score)
        else:
            print(f'ERROR: {res.stderr.decode("utf-8")}')
            sys.exit()
        # cleanup
        os.remove(query)
        # Create output.
        mmdict[mirna[0]] = Mirna(*mirna)

    return mmdict


# Allow boolean argument parsing
def str2bool(v):
    if isinstance(v, bool):
        return v
    if v.lower() in ('yes', 'true', 't', 'y', '1'):
        return True
    elif v.lower() in ('no', 'false', 'f', 'n', '0'):
        return False
    else:
        raise argparse.ArgumentTypeError('Boolean value expected.')


# Main function
def main():
    # Print header
    print('\n'+'#'*57)
    print('###'+' '*51+'###')
    print('###   ncOrtho - ortholog search for non-coding RNAs   ###')
    print('###'+' '*51+'###')
    print('#'*57+'\n')

    ##########################################################################################
    # Command line arguments
    ##########################################################################################
    
    # Parse command-line arguments
    # Define global variables
    parser = argparse.ArgumentParser(
        description='Find orthologs of reference miRNAs in the genome of a query species.'
    )
    parser._action_groups.pop()
    required = parser.add_argument_group('Required Arguments')
    optional = parser.add_argument_group('Optional Arguments')
    # covariance models folder
    required.add_argument(
        '-m', '--models', metavar='<path>', type=str, required=True,
        help='Path to directory containing covariance models (.cm)'
    )
    # mirna data
    required.add_argument(
        '-n', '--ncrna', metavar='<path>', type=str, required=True,
        help='Path to Tab separated file with information about the reference miRNAs'
    )
    # output folder
    required.add_argument(
        '-o', '--output', metavar='<path>', type=str, required=True,
        help='Path to the output directory'
    )
    # query genome
    required.add_argument(
        '-q', '--query', metavar='<.fa>', type=str, required=True,
        help='Path to query genome in FASTA format'
    )
    # reference genome
    required.add_argument(
        '-r', '--reference', metavar='<.fa>', type=str, required=True,
        help='Path to reference genome in FASTA format'
    )
    ##########################################################################
    # Optional Arguments
    ##########################################################################
    # query_name
    optional.add_argument(
        '--queryname', metavar='str', type=str, nargs='?', const='', default='',
        help=(
            'Name for the output directory (RECOMMENDED) '
        )
    )
    # cpu, use maximum number of available cpus unless specified otherwise
    optional.add_argument(
        '--cpu', metavar='int', type=int,
        help='Number of CPU cores to use (Default: all available)', nargs='?',
        const=mp.cpu_count(), default=mp.cpu_count()
    )
    # bit score cutoff for cmsearch hits
    optional.add_argument(
        '--cm_cutoff', metavar='float', type=float,
        help='CMsearch bit score cutoff, given as ratio of the CMsearch bitscore '
             'of the CM against the refernce species (Default: 0.5)', nargs='?', const=0.5, default=0.5
    )
    # length filter to prevent short hits
    optional.add_argument(
        '--minlength', metavar='float', type=float,
        help='CMsearch hit in the query species must have at '
             'least the length of this value times the length of the refernce pre-miRNA (Default: 0.7)',
        nargs='?', const=0.7, default=0.7
    )
    optional.add_argument(
        '--heuristic', type=str2bool, metavar='True/False', nargs='?', const=True, default=True,
        help=(
            'Perform a BLAST search of the reference miRNA in the query genome to identify '
            'candidate regions for the CMsearch. Majorly improves speed. (Default: True)'
        )
    )
    optional.add_argument(
        '--heur_blast_evalue', type=float, metavar='float', nargs='?', const=0.5, default=0.5,
        help=(
            'Evalue filter for the BLASTn search that '
            'determines candidate regions for the CMsearch when running ncOrtho in heuristic mode. '
            '(Default: 0.5) (Set to 10 to turn off)'
        )
    )
    optional.add_argument(
        '--heur_blast_length', type=float, metavar='float', nargs='?', const=0.5, default=0.5,
        help=(
            'Length cutoff for BLASTN search with which candidate regions for the CMsearch are identified.'
            'Cutoff is given as ratio of the reference pre-miRNA length '
            '(Default: 0.5) (Set to 0 to turn off)'
        )
    )
    optional.add_argument(
        '--cleanup', type=str2bool, metavar='True/False', nargs='?', const=True, default=True,
        help=(
            'Cleanup temporary files (Default: True)'
        )
    )
    optional.add_argument(
        '--refblast', type=str, metavar='<path>', nargs='?', const='', default='',
        help=(
            'Path to BLASTdb of the reference species'
        )
    )
    optional.add_argument(
        '--queryblast', type=str, metavar='<path>', nargs='?', const='', default='',
        help=(
            'Path to BLASTdb of the query species'
        )
    )
    optional.add_argument(
        '--maxcmhits', metavar='int', nargs='?', const=None, default=None,
        help=(
            'Maximum number of cmsearch hits to examine. '
            'Decreases runtime significantly if reference miRNA in genomic repeat region. '
            'Set to empty variable to disable (i.e. --maxcmhits=None, default)'
        )
    )
    # use dust filter?
    optional.add_argument(
        '--dust', metavar='yes/no', type=str,
        help='Use BLASTn dust filter during re-BLAST. Greatly decreases runtime if reference miRNA(s) '
             'are located in repeat regions. '
             'However ncOrtho will also not identify orthologs for these miRNAs',
        nargs='?',
        const='no', default='no'
    )
    # check Co-ortholog-ref
    optional.add_argument(
        '--checkCoorthologsRef', type=str2bool, metavar='True/False', nargs='?', const=False, default=False,
        help=(
            'If the re-blast does not identify the original reference miRNA sequence as best hit,'
            'ncOrtho will check whether the best blast '
            'hit is likely a co-ortholog of the reference miRNA relative to the search taxon. '
            'NOTE: Setting this flag will substantially increase'
            'the sensitivity of HaMStR but most likely affect also the specificity, '
            'especially when the search taxon is evolutionarily only very'
            'distantly related to the reference taxon (Default: False)'
        )
    )
    ##########################################################################################
    # Parse arguments
    ##########################################################################################

    # Show help when no arguments are added.
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)
    else:
        args = parser.parse_args()
    
    # Check if computer provides the desired number of cores.
    available_cpu = mp.cpu_count()
    if args.cpu > available_cpu:
        raise ValueError('The provided number of CPU cores is higher than the number available on this system')
    else:
        cpu = args.cpu

    # mandatory
    mirnas = args.ncrna
    models = args.models
    output = args.output
    query = args.query
    reference = args.reference

    # optional
    try:
        max_hits = int(args.maxcmhits)
    except TypeError:
        max_hits = None
    cm_cutoff = args.cm_cutoff
    checkCoorthref = args.checkCoorthologsRef
    cleanup = args.cleanup
    heuristic = (args.heuristic, args.heur_blast_evalue, args.heur_blast_length)
    msl = args.minlength
    refblast = args.refblast
    qblast = args.queryblast
    dust = args.dust.strip()

    ##########################################################################################
    # Starting argument checks
    ##########################################################################################

    # Test if optional query name was given
    if args.queryname:
        qname = args.queryname
    else:
        qname = '.'.join(query.split('/')[-1].split('.')[:-1])
    # make folder for query in output dir
    if not output.split('/')[-1] == qname:
        output = f'{output}/{qname}'

    # Test if input files exist
    all_files = [mirnas, reference, query]
    for pth in all_files:
        if not os.path.isfile(pth):
            raise ValueError(f'{pth} is not a file')
    if not os.path.isdir(models):
        raise ValueError(f'Directory with covariance models does not exist at: {models}')

    # create symbolic links to guarantee writing permissions for pyfaidx index
    q_data = '{}/data'.format(output)
    if not os.path.isdir(q_data):
        os.makedirs(q_data)
    qlink = f'{q_data}/{qname}.fa'
    try:
        os.symlink(query, qlink)
    except FileExistsError:
        pass

    # check if reference BLASTdb was given as input
    if refblast:
        # check if BLASTdb exists
        if not check_blastdb(refblast):
            raise ValueError('# Reference BLASTdb not found at: {}'.format(refblast))
    else:
        # check if refblast exists at location of reference genome
        if not check_blastdb(reference):
            refname = reference.split('/')[-1]
            refblast = f'{q_data}/{refname}'
            # check if BLASTdb exists already. Otherwise create it
            if not check_blastdb(refblast):
                make_blastndb(reference, refblast)
        else:
            refblast = reference
    # check if query BLASTdb was given as input (only used in heuristic mode)
    if qblast and heuristic[0]:
        # check if qblast exists
        if not check_blastdb(qblast):
            raise ValueError('# Query BLASTdb not found at: {}'.format(qblast))
    elif heuristic[0]:
        qblast = qlink
        # check if already exists
        if not check_blastdb(qlink):
            print('# Creating query Database')
            make_blastndb(query, qlink)

    ###############################################################################################
    # Main
    ###############################################################################################

    # Print out query
    print('### Starting ncOrtho run for {}\n'.format(qname))

    # Create miRNA objects from the list of input miRNAs.
    mirna_dict = mirna_maker(mirnas, models, output, msl)

    outpath = '{0}/{1}_orthologs.fa'.format(output, qname)
    log_file = f'{output}/{qname}.log'
    with open(log_file, 'w') as log, open(outpath, 'w') as of:
        log.write(f'# miRNA\tStatus\n')
        # Identify ortholog candidates.
        for mirna in mirna_dict:
            restricted = False
            sys.stdout.flush()
            mirna_data = mirna_dict[mirna]
            # mirna objects for which no CM was found are empty
            if not mirna_data:
                log.write(f'{mirna}\tNo CM\n')
                continue
            print(f'\n### {mirna}')

            # Create output folder, if not existent.
            if not heuristic[0] or not cleanup:
                outdir = '{}/{}'.format(output, mirna)
            else:
                outdir = '{}'.format(output)
            if not os.path.isdir(outdir):
                os.makedirs(outdir)

            # start cmsearch
            cm_results, exitstatus = cmsearcher(
                mirna_data, cm_cutoff, cpu, msl, models, qlink, qblast, outdir, cleanup, heuristic
            )

            # Extract sequences for candidate hits (if any were found).
            if not cm_results:
                print('# {} for {}'.format(exitstatus, mirna))
                log.write(f'{mirna}\t{exitstatus}\n')
                continue
            elif max_hits:
                if len(cm_results) > max_hits:
                    print('# Maximum CMsearch hits reached. Restricting to best {} hits'.format(max_hits))
                    cm_results = {k: cm_results[k] for k in list(cm_results.keys())[:max_hits]}
                    restricted = True
            # get sequences for each CM result
            gp = GenomeParser(qlink, cm_results.values())
            candidates = gp.extract_sequences()
            nr_candidates = len(candidates)
            if nr_candidates == 1:
                print(
                    '# Covariance model search successful, found 1 '
                    'ortholog candidate.'
                )
            else:
                print(
                    '# Covariance model search successful, found {} '
                    'ortholog candidates.'
                    .format(nr_candidates)
                )
            print('# Evaluating candidates.')

            # Perform reverse BLAST test to verify candidates, stored in
            # a list (accepted_hits).
            reblast_hits = {}
            for candidate in candidates:
                sequence = candidates[candidate]
                # print('# Starting reverse blast for {}'.format(candidate))
                blast_command = (
                    'blastn -task blastn -db {0} -num_threads {1} -dust {2} -outfmt "6 qseqid sseqid pident '
                    'length mismatch gapopen qstart qend sstart send evalue bitscore sseq"'
                        .format(refblast, cpu, dust)
                )
                blast_command = (
                    'blastn -task blastn -db {0} -num_threads {1} -dust {2} -outfmt "6 sseqid '
                    'sstart send sstrand bitscore"'
                    .format(refblast, cpu, dust)
                )
                blast_call = sp.Popen(
                    blast_command, shell=True, stdin=sp.PIPE,
                    stdout=sp.PIPE, stderr=sp.PIPE, encoding='utf8'
                )
                # run BLAST
                res, err = blast_call.communicate(sequence)
                if err:
                    print(f'ERROR: {err}')
                    continue
                if not cleanup:
                    blast_output = '{0}/reBLAST_{1}.out'.format(outdir, candidate)
                    with open(blast_output, 'w') as bh:
                        for line in res:
                            bh.write(line)
                # parse BLASTn results
                blast_output = [line.split() for line in res.split('\n')]
                if not blast_output[0]:
                    print('No re-BLAST hits')
                    continue

                # BlastParser will read the temp_fasta file
                bp = BlastParser(mirna_data, blast_output, msl)
                # the parse_blast_output function will return True if the candidate is accepted
                if bp.evaluate_besthit():
                    print('Found best hit')
                    reblast_hits[candidate] = sequence
                elif checkCoorthref:
                    print('Best hit differs from reference sequence! Doing further checks\n')
                    # coorth_out = '{0}/{1}_coorth'.format(outdir, candidate)
                    if bp.check_coortholog_ref(sequence, outdir):
                        reblast_hits[candidate] = sequence
                else:
                    print('Best hit does not overlap with miRNA location')

            # Write output file if at least one candidate got accepted.
            if reblast_hits:
                nr_accepted = len(reblast_hits)
                if nr_accepted == 1:
                    print('# ncOrtho found 1 verified ortholog.')
                    out_dict = reblast_hits
                else:
                    print(
                        '# ncOrtho found {} potential co-orthologs.'
                        .format(nr_accepted)
                    )
                    print('# Evaluating distance between candidates to verify co-orthologs')
                    rbp = ReBlastParser(mirna_data, reblast_hits)
                    out_dict = rbp.verify_coorthologs(outdir)
                    if len(out_dict) == 1:
                        print('ncOrtho found 1 verified ortholog')
                    else:
                        print(
                            '# ncOrtho found {} verified co-orthologs.'
                                .format(len(out_dict))
                        )

                for hit in out_dict:
                    cmres = list(cm_results[hit])
                    cmres.insert(0, qname)
                    cmres = [str(entry) for entry in cmres]
                    header = '|'.join(cmres)
                    of.write('>{0}\n{1}\n'.format(header, out_dict[hit]))

                log.write(f'{mirna}\tSUCESS\n')
            else:
                print(
                    '# None of the candidates for {} could be verified.'
                    .format(mirna)
                )
                if restricted:
                    log.write(f'{mirna}\tNo re-BLAST after restricting to {max_hits} CMsearch hits\n')
                else:
                    log.write(f'{mirna}\tNo re-BLAST\n')
        print('\n### ncOrtho is finished!')


if __name__ == "__main__":
    main()
