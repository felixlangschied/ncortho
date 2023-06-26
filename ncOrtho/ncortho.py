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
import sys
from time import time
from tqdm import tqdm
import logging
from pyfiglet import Figlet
from importlib.metadata import version

# Internal ncOrtho modules
try:
    from blastparser import ReBlastParser
    from genparser import GenomeParser
    from cmsearch import model_search
    from utils import check_blastdb, make_blastndb, write_output, str2bool
    from reblast import perform_reblast
    from rna_object import rna_maker
except ModuleNotFoundError:
    from ncOrtho.blastparser import ReBlastParser
    from ncOrtho.genparser import GenomeParser
    from ncOrtho.cmsearch import model_search
    from ncOrtho.utils import check_blastdb, make_blastndb, write_output, str2bool
    from ncOrtho.reblast import perform_reblast
    from ncOrtho.rna_object import rna_maker


# Main function
def main():
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
        '--sensitive_heuristic', type=float, metavar='True/False', nargs='?', const=False, default=False,
        help=(
            'If no candidate region is found via BLASTn, search with CM in full query genome (Default: False)'
        )
    )
    optional.add_argument(
        '--cleanup', type=str2bool, metavar='True/False', nargs='?', const=True, default=True,
        help=(
            'Cleanup temporary files (Default: True)'
        )
    )
    optional.add_argument(
        '--phmm', type=str2bool, metavar='True/False', nargs='?', const=False, default=False,
        help=(
            'Use pHMM instead of CM for ortholog search (Default: False)'
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
             'However ncOrtho will also not identify orthologs for these miRNAs (Default: no)',
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
        )
    )

    # print header
    custom_fig = Figlet(font='stop')
    print(custom_fig.renderText('ncOrtho')[:-3], flush=True)
    v = version('ncOrtho')
    print(f'Version: {v}\n', flush=True)
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
    output = os.path.realpath(args.output)
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
    heuristic = (args.heuristic, args.heur_blast_evalue, args.heur_blast_length, args.sensitive_heuristic)
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
    if not os.path.isdir(output):
        os.mkdir(output)
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
    if not os.path.islink(qlink):
        os.symlink(query, qlink)

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
            print('# Creating query Database', flush=True)
            make_blastndb(query, qlink)

    ###############################################################################################
    # Main
    ###############################################################################################

    # Print out query
    print('# Starting ncOrtho run for {}\n'.format(qname), flush=True)

    # Create miRNA objects from the list of input miRNAs.
    mirna_dict = rna_maker(mirnas, models, args.phmm, cpu)

    # setup logging
    shortlog = ['# miRNA\tSeconds\tStatus\n']
    logout = f'{output}/{qname}_extended.log'
    if os.path.isfile(logout):
        os.remove(logout)
    logger = logging.getLogger('ncortho')
    logger.setLevel(logging.DEBUG)
    fh = logging.FileHandler(logout)
    fh.setLevel(logging.DEBUG)
    logger.addHandler(fh)
    sucess_count = 0

    print('\n# Ortholog search', flush=True)
    mirna_orthologs = []
    for mirna in tqdm(mirna_dict, file=sys.stdout):
        st = time()
        restricted = False
        sys.stdout.flush()
        mirna_data = mirna_dict[mirna]
        # mirna objects for which no CM was found are empty
        if not mirna_data:
            et = time()
            elapsed_time = str(et - st)
            shortlog.append(f'{mirna}\t{elapsed_time}\tNo CM\n')
            continue
        logger.info(f'\n### {mirna}')

        # Create output folder, if not existent.
        if not heuristic[0] or not cleanup:
            outdir = '{}/{}'.format(output, mirna)
        else:
            outdir = '{}'.format(output)
        if not os.path.isdir(outdir):
            os.makedirs(outdir)

        # start cmsearch
        cm_results, exitstatus = model_search(
            mirna_data, cm_cutoff, cpu, msl, models, qlink, qblast, outdir, cleanup, heuristic, args.phmm
        )

        # Extract sequences for candidate hits (if any were found).
        if not cm_results:
            logger.info('# {} for {}'.format(exitstatus, mirna))
            et = time()
            elapsed_time = str(et - st)
            shortlog.append(f'{mirna}\t{elapsed_time}\t{exitstatus}\n')
            continue
        elif max_hits:
            if len(cm_results) > max_hits:
                logger.warning('# Maximum CMsearch hits reached. Restricting to best {} hits'.format(max_hits))
                cm_results = {k: cm_results[k] for k in list(cm_results.keys())[:max_hits]}
                restricted = True
        # get sequences for each CM result
        gp = GenomeParser(qlink, cm_results.values())
        candidates = gp.extract_sequences()
        # print(candidates)
        logger.info(f'Covariance model search successful, found {len(candidates)} ortholog candidate(s).')

        # Perform reverse BLAST test to verify candidates, stored in
        reblast_hits = {}
        for candidate in candidates:
            sequence = candidates[candidate]
            if list(sequence).count('N') >= 0.9 * len(sequence):
                continue
            # print('# Starting reverse blast for {}'.format(candidate))
            reblasthit = perform_reblast(
                sequence, refblast, cpu, outdir, candidate, mirna_data, dust, msl, cleanup, checkCoorthref
            )
            if reblasthit:
                reblast_hits[candidate] = reblasthit

        # Write output file if no candidate got accepted.
        if not reblast_hits:
            et = time()
            elapsed_time = str(et - st)
            logger.info(f'# None of the candidates for {mirna} could be verified.')
            if restricted:
                shortlog.append(f'{mirna}\t{elapsed_time}\tNo re-BLAST after restricting to {max_hits} CMsearch hits\n')
            else:
                shortlog.append(f'{mirna}\t{elapsed_time}\tNo re-BLAST\n')
            continue

        # Write output file if at least one candidate got accepted.
        nr_accepted = len(reblast_hits)
        if nr_accepted == 1:
            logger.info('# ncOrtho found 1 verified ortholog.')
            out_dict = reblast_hits
        else:
            logger.info(
                '# ncOrtho found {} potential co-orthologs.'
                .format(nr_accepted)
            )
            logger.info('# Evaluating distance between candidates to verify co-orthologs')
            rbp = ReBlastParser(mirna_data, reblast_hits)
            out_dict = rbp.verify_coorthologs(outdir)
            if len(out_dict) == 1:
                logger.info('ncOrtho found 1 verified ortholog')
            else:
                logger.info(f'# ncOrtho found {len(out_dict)} verified co-orthologs.')

        for hit in out_dict:
            cmres = list(cm_results[hit])
            cmres.insert(0, qname)
            cmres = [str(entry) for entry in cmres]
            header = '|'.join(cmres)
            mirna_orthologs.append('>{0}\n{1}\n'.format(header, out_dict[hit]))

        et = time()
        elapsed_time = str(round(et - st, 3))
        shortlog.append(f'{mirna}\t{elapsed_time}\tSUCESS\n')
        sucess_count += 1

    write_output(mirna_orthologs, f'{output}/{qname}_orthologs.fa')
    write_output(shortlog, f'{output}/{qname}_summary.log')
    # logging.basicConfig(filename=f'{output}/{qname}_extended.log', level=logging.DEBUG)

    print(f'\n### ncOrtho found orthologs for {sucess_count} out of {len(mirna_dict)} ncRNAs', flush=True)


if __name__ == "__main__":
    main()

