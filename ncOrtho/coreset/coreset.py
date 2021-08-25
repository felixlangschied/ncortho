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
import pyfaidx
import subprocess as sp
import sys
import yaml


from ncOrtho.coreset.createcm import CmConstructor
from ncOrtho.coreset.core_reblast import blastsearch
from ncOrtho.coreset.coreset_utils import gff_parser
from ncOrtho.coreset.coreset_utils import gtf_parser
from ncOrtho.coreset.coreset_utils import ortho_search

# from createcm import CmConstructor
# from core_reblast import blastsearch
# from coreset_utils import gff_parser
# from coreset_utils import gtf_parser
# from coreset_utils import ortho_search


###############################################################################
def main():
    # Print header
    print('\n' + '#' * 43)
    print('###' + ' ' * 37 + '###')
    print('###   ncOrtho - core set construction   ###')
    print('###' + ' ' * 37 + '###')
    print('#' * 43 + '\n')

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
             'Will not build models from repeat regions. (Default: yes)',
        nargs='?',
        const='yes', default='yes'
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
             'compared to the IDs in the pairwise orthologs file (default:ID) '
             'Options: ID, Name, GeneID, gene_id',
        nargs='?', const='ID=', default='ID'
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
        print(
            '# Error: The provided number of CPU cores is higher than the '
            'number available on this system. Exiting...'
        )
        sys.exit(1)
    else:
        cpu = args.threads

    # parse mandatory arguments
    mirna_path = args.ncrna
    output = args.output

    param_file = args.parameters
    path_dict = {}  # {<name>: {'orthologs': <path>, 'genome': <path>, 'annotation': <path>}}
    with open(param_file, 'r') as param_handle:
        params = yaml.load_all(param_handle, Loader=yaml.FullLoader)
        for entry in params:
            in_type = entry.pop('type')
            if in_type == 'reference':
                ref_annot_path = entry['annotation']
                ref_genome = entry['genome']
            elif in_type == 'core':
                species = entry.pop('name')
                path_dict[species] = entry

    # parse optional arguments
    mgi = args.mgi
    dust = args.dust
    id_t = args.idtype
    create_model = args.create_model

    # check if files exist
    all_paths = []
    for sd in path_dict:
        all_paths.extend(list(path_dict[sd].values()))
    all_paths.extend([ref_annot_path, ref_genome])
    for cp in all_paths:
        if not os.path.isfile(cp):
            print(f'ERROR: {cp} does not exist')
            sys.exit()

    ###############################################################################
    ortho_dict = {}
    neighbor_dict = {}
    mirna_dict = {}

    print('# Reading pairwise orthologs')
    for taxon in path_dict:
        with open(path_dict[taxon]['orthologs']) as omafile:
            oma_lines = omafile.readlines()
            orthologs = {
                ref: core for (ref, core) in [
                    (line.split()[0], line.split()[1])
                    for line in oma_lines if not line.startswith('#')
                ]
            }
            ortho_dict[taxon] = orthologs
    print('Done')

    print('# Reading miRNA data')
    with open(mirna_path) as mirfile:
        mirnas = [
            line.split() for line in mirfile.readlines()
            if not line.startswith('#')
        ]
    print('Done')

    print('# Reading reference annotation')
    ft = ref_annot_path.split('.')[-1]
    if ft == 'gtf':
        ref_dict = gtf_parser(ref_annot_path)
    elif ft in ['gff3', 'gff']:
        ref_dict = gff_parser(ref_annot_path, id_t)
    else:
        print('ERROR: File type "{}" not valid as reference annotation'.format(ft))
        sys.exit()
    print('Done')

    # Determine the position of each miRNA and its neighboring gene(s)
    for mirna in mirnas:
        sys.stdout.flush()
        mirid = mirna[0]
        print('### {0} ###'.format(mirid))
        # Check if output folder exists or create it otherwise
        if not os.path.isdir('{}/{}'.format(output, mirid)):
            sp.call('mkdir -p {}/{}'.format(output, mirid), shell=True)
        ### Workaround for differing naming conventions in miRBase and Ensembl
        if 'chr' in mirna[1]:
            chromo = mirna[1].split('chr')[1]
        else:
            chromo = mirna[1]
        ##
        start = int(mirna[2])
        end = int(mirna[3])
        strand = mirna[4]

        ### find left neighbor or check if located inside gene
        ### chr_dict[contig][i] = (geneid, start, end, strand)
        ###############################################################################
        # case 1): there is no protein-coding gene on the same contig as the miRNA,
        # so there can be no neighbors (should only occur in highly fragmented
        # assemblies)
        if not chromo in ref_dict.keys():
            print(
                'WARNING: No protein-coding genes found on contig "{0}". '
                'Synteny around {1} cannot be established.\n'
                'Make sure that the contig identifiers of the miRNA input file '
                'match the ones in the reference annotation file'
                    .format(chromo, mirid)
            )
            continue

        # case 2): miRNA is located left of the first gene and hence has no left
        # neighbor, the first gene is therefore by default the right neighbor
        if end < int(ref_dict[chromo][1][1]):
            print(
                'There is no left neighbor of {0}, since it is located at the '
                'start of contig {1}.'.format(mirid, chromo)
            )
            print(
                '{0} is the right neighbor of {1}.'
                    .format(ref_dict[chromo][1][0], mirid)
            )
            continue

        # case 3): miRNA is located right to the last gene, so the last gene is the
        # left neighbor and there cannot be a right neighbor
        elif start > int(ref_dict[chromo][len(ref_dict[chromo])][2]):
            print(
                '{0} is the left neighbor of {1}.'
                    .format(ref_dict[chromo][len(ref_dict[chromo])][0], mirid)
            )
            print(
                'There is no right neighbor of {0}, since it is located at the'
                ' end of contig {1}.'.format(mirid, chromo)
            )
            continue

        # case 4): miRNA is located either between two genes or overlapping with (an
        # intron of) a gene, either on the same or the opposite strand
        ###############################################################################
        else:
            solved = False
            for i, gene in enumerate(ref_dict[chromo]):
                gene_data = ref_dict[chromo][gene]
                # case 4.1): miRNA inside gene
                if (
                        start >= gene_data[1]
                        and end <= gene_data[2]
                        and strand == gene_data[3]
                ):
                    solved = True
                    # c+=1
                    print(
                        '{0} is located inside the gene {1}.'
                            .format(mirid, gene_data[0])
                    )
                    ortho_hits = ortho_search(gene_data[0], ortho_dict)
                    for core_tax in ortho_hits:
                        try:
                            neighbor_dict[core_tax][mirid] = (
                                ('inside', ortho_hits[core_tax])
                            )
                        except:
                            neighbor_dict[core_tax] = (
                                {mirid: ('inside', ortho_hits[core_tax])}
                            )
                    break
                # case 4.2): miRNA opposite of gene
                elif (
                        start >= gene_data[1]
                        and end <= gene_data[2]
                        and strand != gene_data[3]
                ):
                    solved = True
                    print(
                        '{0} is located opposite of the gene {1}.'
                            .format(mirid, gene_data[0])
                    )
                    ortho_hits = ortho_search(gene_data[0], ortho_dict)
                    for core_tax in ortho_hits:
                        try:
                            neighbor_dict[core_tax][mirid] = (
                                ('opposite', ortho_hits[core_tax])
                            )
                        except:
                            neighbor_dict[core_tax] = (
                                {mirid: ('opposite', ortho_hits[core_tax])}
                            )
                    break
                # case 4.3): miRNA between genes
                elif (
                        int(ref_dict[chromo][gene][2]) < start
                        and ref_dict[chromo][gene + 1][1] > end
                ):
                    solved = True
                    ###############################################################################
                    print(
                        '{1} is the left neighbor of {2}.'
                        .format(gene, ref_dict[chromo][gene][0], mirid)
                    )
                    print(
                        '{1} is the right neighbor of {2}.'
                        .format(gene, ref_dict[chromo][gene + 1][0], mirid)
                    )
                    left_hits = ortho_search(gene_data[0], ortho_dict)
                    right_hits = (
                        ortho_search(ref_dict[chromo][gene + 1][0], ortho_dict)
                    )
                    # save only the hits where both genes have orthologs in a species
                    if left_hits:
                        # print(left_hits)
                        # print(right_hits)
                        for taxon in left_hits:
                            if taxon in right_hits:
                                try:
                                    neighbor_dict[taxon][mirid] = (
                                        (
                                            'in-between',
                                            [left_hits[taxon],
                                             right_hits[taxon]]
                                        )
                                    )
                                except:
                                    neighbor_dict[taxon] = (
                                        {mirid: (
                                            'in-between',
                                            [left_hits[taxon],
                                             right_hits[taxon]]
                                        )}
                                    )
                            else:
                                print('Orthologs were not found for both flanking genes')
                    break
            if not solved:
                print('Unable to resolve synteny for {}.'.format(mirid))

    #################################################################################################

    # Search for the coordinates of the orthologs and extract the sequences
    for taxon in neighbor_dict:
        print('\n### Starting synteny analysis for {}'.format(taxon))
        print('# Trying to parse annotation file for {}.'.format(taxon))
        fe = path_dict[taxon]['annotation'].split('.')[-1]
        if fe == 'gtf':
            core_dict = gtf_parser(path_dict[taxon]['annotation'])
        elif fe in ['gff3', 'gff']:
            core_dict = gff_parser(path_dict[taxon]['annotation'], id_t)
        else:
            print('ERROR: File type "{}" not valid as reference annotation'.format(ft))
            sys.exit()

        print('# Loading genome file')
        fasta_path = path_dict[taxon]['genome']
        core_gen_dict = '{}/core_genomes'.format(output)
        if not os.path.isdir(core_gen_dict):
            os.mkdir(core_gen_dict)
        slink = '{}/slink_to_{}'.format(core_gen_dict, taxon)
        try:
            os.symlink(fasta_path, slink)
        except FileExistsError:
            pass
        genome = pyfaidx.Fasta(slink)
        print('Done')

        # print('YOU MADE IT THIS FAR.')
        for mirna in neighbor_dict[taxon]:
            # print(mirna)
            style = neighbor_dict[taxon][mirna][0]
            if style == 'inside' or style == 'opposite':
                print(f'# {mirna} location: {style} of gene')
                ###############################################################################
                try:
                    ortho_data = (
                        core_dict[neighbor_dict[taxon][mirna][1]]
                    )
                except KeyError as e:
                    print('{} not found in annotation file.'.format(e))
                    continue
                positions = list(
                    core_dict[ortho_data[0]][ortho_data[1]][1:4]
                )
                coordinates = [ortho_data[0]] + positions
                seq = (
                    genome[coordinates[0]]
                    [coordinates[1] - 1:coordinates[2]].seq
                )
                # print(seq[0:10])
                try:
                    mirna_dict[mirna][taxon] = seq
                except KeyError:
                    mirna_dict[mirna] = {taxon: seq}
            elif style == 'in-between':
                print(f'# {mirna} location: between genes')
                try:
                    left_data = (
                        core_dict[neighbor_dict[taxon][mirna][1][0]]
                    )
                    right_data = (
                        core_dict[neighbor_dict[taxon][mirna][1][1]]
                    )
                except KeyError as e:
                    print('{} not found in annotation file.'.format(e))
                    continue
                # Test to see if the two orthologs are themselves neighbors where their
                # distance cannot be larger than the selected mgi value. This accounts
                # for insertions in the core species.
                # TODO: Apply mgi also to the reference species to account for insertions
                # in the reference.
                #                 print(f'left_data: {left_data}')
                #                 print(f'right_data: {right_data}')
                if (
                        left_data[0] == right_data[0]
                        and abs(left_data[1] - right_data[1]) <= mgi
                ):
                    # Determine which sequence to include for the synteny-based ortholog search
                    # depending on the order of orthologs. The order of the orthologs in the core
                    # species might be inverted compared to that in the reference species.
                    ###############################################################################
                    if left_data[1] < right_data[1]:
                        # print('left')
                        # print(core_gtf_dict[left_data[0]][left_data[1]])
                        # print(core_gtf_dict[right_data[0]][right_data[1]])
                        contig = left_data[0]
                        # print(contig)
                        seq_start = int(
                            core_dict[left_data[0]][left_data[1]][2]
                        )
                        # print(seq_start)
                        seq_end = (
                            core_dict[right_data[0]][right_data[1]][1]
                        )
                        # print(seq_end)
                        seq = genome[contig][seq_start - 1:seq_end].seq
                        try:
                            mirna_dict[mirna][taxon] = seq
                        except:
                            mirna_dict[mirna] = {taxon: seq}
                    elif right_data[1] < left_data[1]:
                        # print('right')
                        # print(core_gtf_dict[left_data[0]][left_data[1]])
                        # print(core_gtf_dict[right_data[0]][right_data[1]])
                        contig = left_data[0]
                        # print(contig)
                        # seq_start = int(
                        #     core_gtf_dict[left_data[0]][left_data[1]][2]
                        # )
                        # print(seq_start)
                        # seq_end = (
                        #     core_gtf_dict[right_data[0]][right_data[1]][1]
                        # )
                        seq_start = int(
                            core_dict[right_data[0]][right_data[1]][2]
                        )
                        # print(seq_start)
                        seq_end = (
                            core_dict[left_data[0]][left_data[1]][1]
                        )
                        # print(seq_end)
                        seq = genome[contig][seq_start - 1:seq_end].seq
                        try:
                            mirna_dict[mirna][taxon] = seq
                        except:
                            mirna_dict[mirna] = {taxon: seq}
                    print('Synteny fulfilled.')
                else:
                    print(
                        'No shared synteny for {} in {}.'
                            .format(mirna, taxon)
                    )
                    continue
                    # print(left_data)
                    # print(right_data)
            else:
                print('## Neither inside, opposite, nor in-between')
                # print(neighbor_dict[taxon][mirna])
                continue
            print('Candidate region found')

    # Write output file
    for mirna in mirna_dict:
        # print('{0}/{1}/{1}.fa'.format(output, mirna))
        with open('{0}/{1}/{1}.fa'.format(output, mirna), 'w') as outfile:
            for core_taxon in mirna_dict[mirna]:
                outfile.write(
                    '>{0}\n{1}\n'
                    .format(core_taxon, mirna_dict[mirna][core_taxon])
                )
    if not mirna_dict:
        print('\nWARNING: No syntenic regions found in the core species')

    for mirna in mirnas:
        print('\n### Starting reciprocal BLAST process')
        sto_path = blastsearch(mirna, ref_genome, output, cpu, dust)

        if create_model == 'yes' and sto_path is not None:
            print('### Starting to construct covariance model from alignment')
            model_dir = f'{output}/CMs'
            if not os.path.isdir(model_dir):
                os.mkdir(model_dir)
            model_out = f'{model_dir}/{mirna[0]}.cm'
            if not os.path.isfile(model_out):
                # Initiate covariance model construction and calibration.
                cmc = CmConstructor(sto_path, model_dir, mirna[0], cpu)
                # Construct the model.
                cmc.construct()
                # Calibrate the model.
                cmc.calibrate()
            else:
                print(f'Model of {mirna[0]} already found at {model_dir}. Nothing done..')

    print('\n### Construction of core set finished')


if __name__ == '__main__':
    main()
