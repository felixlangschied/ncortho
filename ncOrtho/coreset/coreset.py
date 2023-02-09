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
import shutil

try:
    from createcm import CmConstructor
    from core_reblast import blastsearch
    from locate_position import categorize_mirna_position
    from synteny import analyze_synteny
    import coreset_utils as cu
except ModuleNotFoundError:
    from ncOrtho.coreset.createcm import CmConstructor
    from ncOrtho.coreset.core_reblast import blastsearch
    from ncOrtho.coreset.synteny import analyze_synteny
    from ncOrtho.coreset.locate_position import categorize_mirna_position
    import ncOrtho.coreset.coreset_utils as cu


###############################################################################
def main():
    # Print header
    print('\n' + '#' * 43, flush=True)
    print('###' + ' ' * 37 + '###', flush=True)
    print('###   ncOrtho - core set construction   ###', flush=True)
    print('###' + ' ' * 37 + '###', flush=True)
    print('#' * 43 + '\n', flush=True)

    # Parse command-line arguments
    # Define global variables
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
    # cpu, use maximum number of available cpus if not specified otherwise
    optional.add_argument(
        '--threads', metavar='int', type=int,
        help='Number of CPU cores to use (Default: All available)', nargs='?',
        const=mp.cpu_count(), default=mp.cpu_count()
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
    #
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
    output = args.output

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
        os.makedirs(tmpout)

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
            print(f'Warning: Could not localize {mirid}')
            continue
        mirna_positions[mirid] = (syntype, core_orthos)

    print('### Identifying syntenic regions in core species', flush=True)
    syntenyregion_per_mirna = analyze_synteny(core_dict, mirna_positions, tmpout, args.idtype, args.mgi, verbose)

    if not syntenyregion_per_mirna:
        raise ValueError('No regions of conserved synteny found in any core species')

    for mirid, fastalist in syntenyregion_per_mirna.items():
        if not fastalist:
            print(f'Warning: No syntenic region found in any core species for {mirid}. '
                  f'Make sure that the IDs in the annotation file match the ones in the orthologs file')
            continue
        miroutdir = f'{tmpout}/{mirid}'
        if not os.path.isdir(miroutdir):
            os.mkdir(miroutdir)
        with open(f'{miroutdir}/synteny_regions_{mirid}.fa', 'w') as of:
            for line in fastalist:
                of.write(line)

    #################################################################################################
    print('\n### Starting Ortholog search')
    for mirna in mirnas:
        mirid = mirna[0]
        sto_path = blastsearch(mirna, ref_paths['genome'], tmpout, cpu, dust, verbose)
        if create_model == 'yes' and sto_path is not None:
            print(f'# Starting to construct covariance model for {mirid}', flush=True)
            model_out = f'{output}/{mirid}.cm'
            if not os.path.isfile(model_out):
                # Initiate covariance model construction and calibration.
                cmc = CmConstructor(sto_path, output, mirid, cpu)
                # Construct the model.
                cmc.construct()
                # Calibrate the model.
                cmc.calibrate()
            else:
                print(f'Model of {mirid} already found at {output}. Nothing done..', flush=True)
    print('\n### Construction of core set finished')


if __name__ == '__main__':
    main()
