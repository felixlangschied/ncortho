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

# Create a core set of orthologs
# Find the corresponding syntenic regions in reference and core species
# Search for core orthologs by reciprocal BLAST search
# Create Stockholm structural alignment

# required:
# reference microRNA data (sequence, coordinates)
# reference taxon: genome, blastdb, gtf file with gene coordinates
# core set taxa: genome, gtf file, pairwise orthologs

import argparse
import multiprocessing as mp
import os
import sys
from pyfiglet import Figlet
from importlib.metadata import version
import shutil

try:
    from createcm import create_cm, create_phmm
    from core_reblast import blastsearch
    from locate_position import categorize_mirna_position
    from synteny import analyze_synteny
    import coreset_utils as cu
    from core_mmseq import coreorthologs_from_mmseq
except ModuleNotFoundError:
    from ncOrtho.coreset.createcm import create_cm, create_phmm
    from ncOrtho.coreset.core_reblast import blastsearch
    from ncOrtho.coreset.locate_position import categorize_mirna_position
    from ncOrtho.coreset.synteny import analyze_synteny
    import ncOrtho.coreset.coreset_utils as cu
    from ncOrtho.coreset.core_mmseq import coreorthologs_from_mmseq


###############################################################################
def main():
    parser = argparse.ArgumentParser(
        description=(
            'Build Covariance models of reference miRNAs from core set of orthologs.'
        )
    )
    parser._action_groups.pop()
    required = parser.add_argument_group('Required Arguments')
    optional = parser.add_argument_group('Optional Arguments')
    # input file path
    required.add_argument('-p', '--parameters', metavar='<path>', type=str, required=True,
                          help='Path to the parameters file in yaml format')
    # mirna data
    required.add_argument(
        '-n', '--ncrna', metavar='<path>', type=str, required=True,
        help='Path to tab separated file of reference miRNAs information '
    )
    # output folder
    required.add_argument(
        '-o', '--output', metavar='<path>', type=str, required=True,
        help='Path for the output folder'
    )
    # OPTIONAL VARIABLES
    # cpu
    optional.add_argument(
        '--threads', metavar='int', type=int,
        help='Number of CPU cores to use (Default: 4)', nargs='?',
        const=4, default=4
    )
    # Maximum gene insertions
    optional.add_argument(
        '--mgi', metavar='int', type=int,
        help='Maximum number of gene insertions in the core species (Default: 3)', nargs='?',
        const=3, default=3
    )
    # use dust filter?
    optional.add_argument(
        '--dust', metavar='yes/no', type=str,
        help='Use BLASTn dust filter. '
             'Will not build models from repeat regions if "yes". (Default: no)',
        nargs='?',
        const='no', default='no'
    )
    optional.add_argument(
        '--create_model', metavar='yes/no', type=str,
        help='set to "no" if you only want to create the alignment. (Default: yes)', nargs='?',
        const='yes', default='yes'
    )
    optional.add_argument(
        '--phmm', metavar='yes/no', type=str,
        help='set to "yes" if you want to create a pHMM instead of a CM (Default: no)', nargs='?',
        const='no', default='no'
    )
    optional.add_argument(
        '--rcoffee', metavar='yes/no', type=str,
        help='set to "no" to use default "t_coffee" instead of "r_coffee" (Default: yes)', nargs='?',
        const='yes', default='yes'
    )
    optional.add_argument(
        '--redo', metavar='yes/no', type=str,
        help='set to "no" to reuse models found at output destination (Default: yes)', nargs='?',
        const='yes', default='yes'
    )
    optional.add_argument(
        '--idtype', metavar='str', type=str,
        help='Choose the ID in the reference gff file that is '
             'compared to the IDs in the pairwise orthologs file (default:GeneID) '
             'Options: ID, Name, GeneID, gene_id, CDS',
        nargs='?', const='ID=', default='ID'
    )
    optional.add_argument(
        '--max_anchor_dist', metavar='int', type=int,
        help='Number of additional genes to the left and right '
             'of the reference miRNA that are to be considered as syntenic anchors. (Default: 3)',
        nargs='?', const=3, default=3
    )
    optional.add_argument(
        '--verbose', metavar='yes/no', type=str,
        help='Print additional information (Default: no)',
        nargs='?', const='no', default='no'
    )
    optional.add_argument(
        '--search', metavar='str', type=str,
        help='Sequence similarity search program. BLASTn or MMseq2 (Default: BLASTn, MMseq2 is currently not recommended)',
        nargs='?', const='BLASTn', default='BLASTn'
    )
    optional.add_argument(
        '--refdb', type=str, metavar='<path>', nargs='?', const='', default='',
        help=(
            'Path to BLAST or MMseq2 database of the reference species (recommended when running ncCreate in parallel).'
            'If using MMseq2, consider to create the database and its index on a local disc (e.g. /tmp/) '
            'prior to running ncCreate'
        )
    )
    optional.add_argument(
        '--reblastmode', type=str, metavar='str', nargs='?', const='blastn', default='blastn',
        help=(
            'BLASTn mode employed for establishing the reciprocal best hit '
            'criterion of ortholog candidates in the reference genome. '
            'Choose between "blastn" (slow and sensitive) and "megablast" (faster but less sensitive) (Default: blastn)'
        )
    )

    # print header
    custom_fig = Figlet(font='stop')
    print(custom_fig.renderText('ncOrtho')[:-3], flush=True)
    v = version('ncOrtho')
    print(f'Version: {v}\n', flush=True)
    ###############################################################################

    # Show help when no arguments are added.
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)
    else:
        args = parser.parse_args()

    # Check if computer provides the desired number of cores.
    available_cpu = mp.cpu_count()
    if args.threads > available_cpu:
        raise ValueError('The provided number of CPU cores is higher than the number available on this system')
    else:
        cpu = args.threads

    # parse mandatory arguments
    output = os.path.realpath(args.output)
    if not os.path.isdir(output):
        os.mkdir(output)

    # parse optional arguments
    # mgi = args.mgi
    dust = args.dust
    create_model = args.create_model
    if args.verbose == 'yes':
        verbose = True
    elif args.verbose == 'no':
        verbose = False
    else:
        raise ValueError(f'Unknown value for "verbose": {args.verbose}')
    add_pos_orthos = args.max_anchor_dist

    if args.reblastmode not in ['blastn', 'megablast']:
        raise ValueError(f'Unknown BLASTn modus "{args.reblastmode}". Choose "blastn" or "megablast"')

    # parameters
    core_dict, ref_paths, all_paths = cu.parse_yaml(args.parameters)
    # check if files exist
    for cp in all_paths:
        if not os.path.isfile(cp):
            raise ValueError(f'{cp} does not exist')

    ###############################################################################
    neighbor_dict = {}
    mirna_dict = {}

    # create directory for intermediate files
    tmpout = f'{output}/tmp'
    if not os.path.isdir(tmpout):
        os.mkdir(tmpout)

    print('# Reading pairwise orthologs', flush=True)
    if args.idtype == 'CDS':
        ortho_dict = cu.pairwise_orthologs_from_cds(core_dict, ref_paths['annotation'])
    else:
        ortho_dict = cu.read_pairwise_orthologs(core_dict)

    print('# Reading miRNA data', flush=True)
    mirnas = cu.read_mirnas(args.ncrna)

    print('# Reading reference annotation', flush=True)
    ref_dict = cu.parse_annotation(ref_paths['annotation'], args.idtype)

    print('# Determining the position of each miRNA and its neighboring gene(s)', flush=True)
    mirna_positions = {}
    for mirna in mirnas:
        mirid, chromo, start, end, strand = cu.mirna_position(mirna)
        syntype, core_orthos = categorize_mirna_position(
            mirid, chromo, start, end, strand, ref_dict, ortho_dict, add_pos_orthos, verbose
        )
        if not syntype:
            print(f'Warning: Could not localize {mirid}', flush=True)
            continue
        mirna_positions[mirid] = (syntype, core_orthos)

    print('\n### Identifying syntenic regions in core species', flush=True)
    syntenyregion_per_mirna = analyze_synteny(core_dict, mirna_positions, tmpout, args.idtype, args.mgi, verbose)

    if not syntenyregion_per_mirna:
        raise ValueError('No regions of conserved synteny found in any core species')

    for mirid, fastalist in syntenyregion_per_mirna.items():
        if not fastalist:
            print(f'Warning: No syntenic region found in any core species for {mirid}. '
                  f'Make sure that the IDs in the annotation file match the ones in the orthologs file', flush=True)
            continue
        miroutdir = f'{tmpout}/{mirid}'
        if not os.path.isdir(miroutdir):
            os.mkdir(miroutdir)
        with open(f'{miroutdir}/synteny_regions_{mirid}.fa', 'w') as of:
            for line in fastalist:
                of.write(line)
    del syntenyregion_per_mirna  # free memory

    #################################################################################################
    print('\n### Reciprocal best hit search to find core orthologs and training models', flush=True)
    for mirna in mirnas:
        mirid = mirna[0]
        miroutdir = f'{tmpout}/{mirid}'

        if args.search == 'BLASTn':
            corefile = blastsearch(
                mirna, args.refdb, ref_paths['genome'], tmpout, cpu, dust, args.reblastmode, verbose
            )
        elif args.search == 'MMseq2':
            corefile = coreorthologs_from_mmseq(mirna, args.refdb, ref_paths['genome'], tmpout, cpu, dust, verbose)
        else:
            raise ValueError(f'Unknown search program "{args.search}". Choose between "BLASTn" and "MMseq2"')

        if create_model == 'no':
            shutil.copy(corefile, output)
            continue

        sto_path = cu.make_alignment(miroutdir, mirid, cpu, corefile, args.rcoffee)
        if args.phmm == 'no':
            print(f'# Starting to construct covariance model for {mirid}', flush=True)
            model_out = f'{output}/{mirid}.cm'
            if args.redo == 'no' and not os.path.isfile(model_out):
                print(f'Model of {mirid} already found at {output}. Nothing done..', flush=True)
                continue
            create_cm(sto_path, output, mirid, cpu)

        elif args.phmm == 'yes':
            print(f'# Starting to construct pHMM for {mirid}', flush=True)
            create_phmm(sto_path, output, mirid, cpu)
        else:
            raise ValueError(f'Unknown value "{args.phmm}" for --phmm')

    print('\n### Construction of core set finished', flush=True)


if __name__ == '__main__':
    main()
