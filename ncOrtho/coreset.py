# Create a core set of orthologs
# Find the corresponding syntenic regions in reference and core species
# Search for core orthologs by reciprocal BLAST search
# Create Stockholm structural alignment

#required:
#reference microRNA data (sequence, coordinates)
#reference taxon: genome, blastdb, gtf file with gene coordinates
#core set taxa: genome, gtf file, pairwise orthologs

import argparse
import glob
import multiprocessing as mp
import os
import pyfaidx
import subprocess as sp
import sys
import re

from createcm import CmConstructor
from core_reblast import blastsearch

###############################################################################

# Parse an Ensembl GTF file to store the coordinates for each protein-coding gene in a
# dictionary
def gtf_parser(gtf):
    species_name = gtf.split('/')[-1].split('.')[0]
    chr_dict = {}
    chromo = ''

    #with open(inpath) as infile, open(outpath, 'wb') as outfile:
    with open(gtf) as infile:
        for line in infile:
            if (
                not line.startswith('#')
                and line.split()[2] == 'gene'
                and line.split('gene_biotype')[1].split('\"')[1]
                == 'protein_coding'
            ):
                linedata = line.strip().split('\t')
                contig = linedata[0]
                geneid = linedata[-1].split('\"')[1]
                start = int(linedata[3])
                end = int(linedata[4])
                strand = linedata[6]
                if contig != chromo:
                    i = 1
                    chromo = contig
                chr_dict[geneid] = (contig, i)
                try:
                    chr_dict[contig][i] = (geneid, start, end, strand)
                except:
                    chr_dict[contig] = {i: (geneid, start, end, strand)}
                i += 1
    return chr_dict


def gff_parser(gff, id_type):
    """

    Parameters
    ----------
    gff :   Path to a RefSeq gff3 file
    id  :   type of identifier in the gff3 file that is to be compared against the IDs in the pairwise ortholgs file

    Returns
    -------
    Dictionary containing the location of each gene on a chromosome
    {'gene_id': ('chromosome': position), 'chromosome': {position: ('gene_id', start, end, strand)} }

    """
    chr_dict = {}
    chromo = ''

    # with open(inpath) as infile, open(outpath, 'wb') as outfile:
    with open(gff) as infile:
        for line in infile:
            if (
                    not line.startswith('#')
                    and line.split()[2] == 'gene'
                    and line.split('gene_biotype=')[1].split(';')[0]
                    == 'protein_coding'
            ):
                linedata = line.strip().split('\t')
                if id_type in ['ID', 'Name', 'gene_id']:
                    geneid = linedata[-1].split(f'{id_type}=')[1].split(';')[0]
                elif id_type == 'GeneID':
                    geneid = re.split('[;,]', linedata[-1].split(f'{id_type}:')[1])[0]
                else:
                    print('ERROR: Unknown identifier type "{}"'.format(id_type))
                    sys.exit()

                contig = linedata[0]
                # geneid = linedata[-1].split(id)[1].split(';')[0]
                start = int(linedata[3])
                end = int(linedata[4])
                strand = linedata[6]
                if contig != chromo:
                    i = 1
                    chromo = contig
                chr_dict[geneid] = (contig, i)
                try:
                    chr_dict[contig][i] = (geneid, start, end, strand)
                except:
                    chr_dict[contig] = {i: (geneid, start, end, strand)}
                i += 1
    return chr_dict

###############################################################################
#Try to find the ortholog for a given reference gene in a core set species
def ortho_search(r_gene, ortho_dict):
    orthologs = {}
    for core_taxon in ortho_dict.keys():
        try:
            ortholog = ortho_dict[core_taxon][r_gene]
            orthologs[core_taxon] = ortholog
            print(
                '{0} is the ortholog for {1} in {2}.'
                .format(ortholog, r_gene, core_taxon)
            )
        except:
            print(
                'No ortholog found for {0} in {1}.'
                .format(r_gene, core_taxon)
            )
    return orthologs



###############################################################################
def main():
    #c = 0
    ortho_dict = {}
    mirna_dict = {}
    neighbor_dict = {}

    # Print header
    print('\n'+'#'*43)
    print('###'+' '*37+'###')
    print('###   ncOrtho - core set construction   ###')
    print('###'+' '*37+'###')
    print('#'*43+'\n')

# required arguments
# ref gtf, core gtf, core genomes, oma orthologs, mirnas, output, (cpu), (mip)
#python {} -r {} -c {} -g {} -m {} -o {}    
    # Parse command-line arguments
    # Define global variables
    parser = argparse.ArgumentParser(
        prog='python coreset.py', description='core set construction'
    )
    # mirna data
    parser.add_argument(
        '-n', '--ncrna', metavar='<path>', type=str,
        help='path to your reference micrornas'
    )
    # output folder
    parser.add_argument(
        '-o', '--output', metavar='<path>', type=str,
        help='path for the output folder'
    )
    # reference gtf
    parser.add_argument(
        '-r', '--reference', metavar='<.fa>', type=str,
        help='path to reference GTF file'
    )
    # reference genome
    parser.add_argument(
        '-g', '--genome', metavar='<.fa>', type=str,
        help='path to reference genome'
    )
    # core taxa GTF
    parser.add_argument(
        '-c', '--core', metavar='<path>', type=str,
        help='path to core GTF files'
    )
    # core taxa genome
    parser.add_argument(
        '-q', '--query', metavar='<path>', type=str,
        help='path to core genomes'
    )
    # pairwise orthologs folder
    parser.add_argument(
        '-p', '--pairwise', metavar='<path>', type=str,
        help='path to pairwise orthologs'
    )

    # OPTIONAL VARIABLES
    # cpu, use maximum number of available cpus if not specified otherwise
    parser.add_argument(
        '-t', '--threads', metavar='int', type=int,
        help='number of CPU cores to use', nargs='?',
        const=mp.cpu_count(), default=mp.cpu_count()
    )
    # Maximum gene insertions
    parser.add_argument(
        '-m', '--mgi', metavar='int', type=int,
        help='maximum number of gene insertions', nargs='?',
        const=3, default=3
    )
    # use dust filter?
    parser.add_argument(
        '--dust', metavar='yes/no', type=str,
        help='Use BLASTn dust filter. Decreases number of models created but improves runtime and possibly specificity', nargs='?',
        const='yes', default='yes'
    )
    parser.add_argument(
        '--create_model', metavar='yes/no', type=str,
        help='set to "no" if you only want to create the alignment', nargs='?',
        const='yes', default='yes'
    )
    #
    parser.add_argument(
        '--id_type', metavar='str', type=str,
        help='Choose the ID in the reference gff file that is '
             'compared to the IDs in the pairwise orthologs file. Default:"ID"'
             'Options: "Name", "GeneID", gene_id',
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

    ###TODO: include checks for validity of arguments
    #os.getcwd()
    #os.chdir(path)
    #os.path.exists(path)
    #os.path.isfile(path)
    #os.path.isdir(path)

    mirna_path = args.ncrna
    output = args.output
    query = args.query
    ref_annot_path = args.reference
    ref_genome = args.genome
    oma_paths = glob.glob('{}/*'.format(args.pairwise))
    core_gtf_paths = args.core
    core_fa_paths = args.query
    mgi = args.mgi
    dust = args.dust
    id_t = args.id_type
    create_model = args.create_model

###############################################################################

# Parse the pairwise orthologs
    for oma_path in oma_paths:
        taxon = oma_path.split('/')[-1].split('.')[0]
        with open(oma_path) as omafile:
            oma_lines = omafile.readlines()
            orthologs = {
                ref: core for (ref, core) in [
                    (line.split()[0], line.split()[1])
                    for line in oma_lines
                ]
            }
            ortho_dict[taxon] = orthologs

#Read in the miRNA data
    with open(mirna_path) as mirfile:
        mirnas = [
            line.split() for line in mirfile.readlines()
            if not line.startswith('#')
        ]
    # Read reference annotation
    ft = ref_annot_path.split('.')[-1]
    if ft == 'gtf':
        ref_dict = gtf_parser(ref_annot_path)
    elif ft in ['gff3', 'gff']:
        ref_dict = gff_parser(ref_annot_path, id_t)
    else:
        print('ERROR: File type "{}" not valid as reference annotation'.format(ft))
        sys.exit()


#Determine the position of each miRNA and its neighboring gene(s)
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
###
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
                'There are no protein-coding genes on contig {0}. '
                'Synteny around {1} cannot be established.'
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
            ### case 4.1): miRNA inside gene
                if (
                    start >= gene_data[1]
                    and end <= gene_data[2]
                    and strand == gene_data[3]
                ):
                    solved = True
                    #c+=1
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
            ### case 4.2): miRNA opposite of gene
                elif (
                    start >= gene_data[1]
                    and end <= gene_data[2]
                    and strand != gene_data[3]
                ):
                    solved = True
                    #c+=1
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
            ### case 4.3): miRNA between genes
                elif (
                    int(ref_dict[chromo][gene][2]) < start
                    and ref_dict[chromo][gene+1][1] > end
                ):
                    solved = True
###############################################################################
                    print(
                        '{1} is the left neighbor of {2}.'
                        .format(gene, ref_dict[chromo][gene][0], mirid)
                    )   
                    print(
                        '{1} is the right neighbor of {2}.'
                        .format(gene, ref_dict[chromo][gene+1][0], mirid)
                    )
                    left_hits = ortho_search(gene_data[0], ortho_dict)
                    right_hits = (
                       ortho_search(ref_dict[chromo][gene+1][0], ortho_dict)
                    )
                #save only the hits where both genes have orthologs in a species
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

    # print(neighbor_dict)

### Search for the coordinates of the orthologs and extract the sequences
    # print('# starting now with coordinate search\n')
    for taxon in neighbor_dict:
        print('\n### Starting synteny analysis for {}'.format(taxon))
        print('# Trying to parse annotation file for {}.'.format(taxon))
        gtf_path_list = glob.glob('{0}/*{1}*.gtf'.format(core_gtf_paths, taxon))
        if len(gtf_path_list) > 1:
            print('ERROR: Ambiguous core taxon annotation for {}'.format(taxon))
            sys.exit()
        elif len(gtf_path_list) == 1:
            gtf_path = gtf_path_list[0]
            core_gtf_dict = gtf_parser(gtf_path)
            # print(core_gtf_dict)
        elif len(gtf_path_list) == 0:
            gff_path_list = glob.glob('{0}/*{1}*.gff*'.format(core_gtf_paths, taxon))
            if len(gff_path_list) == 1:
                core_dict = gff_parser(gff_path_list[0], id_t)
            else:
                print('ERROR: Ambiguous core taxon annotation for {}'.format(taxon))
                sys.exit()
        print('Done')
        print('# Loading genome file')
        fasta_path = glob.glob('{0}/{1}*.fa'.format(core_fa_paths, taxon))
        #print(gtf_path)
        #print(fasta_path)
        if len(fasta_path) != 1:
            print('Unable to identify genome file for {}'.format(taxon))
            sys.exit()
        genome = pyfaidx.Fasta(fasta_path[0])
        print('Done')
        # try:
        #     core_gtf_dict = gtf_parser(gtf_path)
        #     print('Done')
        # except:
        #     print('No GTF file found for {}'.format(taxon))
        #     sys.exit()
        # print('YOU MADE IT THIS FAR.')
        for mirna in neighbor_dict[taxon]:
            # print(mirna)
            style = neighbor_dict[taxon][mirna][0]
            if style == 'inside' or style == 'opposite':
                print(f'# {mirna} location: {style} of gene')
###############################################################################
                try:
                    ortho_data = (
                        core_gtf_dict[neighbor_dict[taxon][mirna][1]]
                    )
                except KeyError as e:
                    print('{} not found in annotation file.'.format(e))
                    continue
                positions = list(
                    core_gtf_dict[ortho_data[0]][ortho_data[1]][1:4]
                )
                coordinates = [ortho_data[0]] + positions
                seq = (
                    genome[coordinates[0]]
                    [coordinates[1]-1:coordinates[2]].seq
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
                        core_gtf_dict[neighbor_dict[taxon][mirna][1][0]]
                    )
                    right_data = (
                        core_gtf_dict[neighbor_dict[taxon][mirna][1][1]]
                    )
                except KeyError as e:
                    print('{} not found in annotation file.'.format(e))
                    continue
                # print('#########################')
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
                            core_gtf_dict[left_data[0]][left_data[1]][2]
                        )
                        # print(seq_start)
                        seq_end = (
                            core_gtf_dict[right_data[0]][right_data[1]][1]
                        )
                        # print(seq_end)
                        seq = genome[contig][seq_start-1:seq_end].seq
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
                            core_gtf_dict[right_data[0]][right_data[1]][2]
                        )
                        # print(seq_start)
                        seq_end = (
                            core_gtf_dict[left_data[0]][left_data[1]][1]
                        )
                        # print(seq_end)
                        seq = genome[contig][seq_start-1:seq_end].seq
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

    print('\n### Starting reciprocal BLAST process')
    sto_path = blastsearch(mirna_path, ref_genome, output, cpu, dust)

    if create_model == 'yes':
        print('### Starting to construct covariance model from alignment')
        model_out = f'{output}/CMs'
        if not os.path.isdir(model_out):
            os.mkdir(model_out)
        name = sto_path.split('/')[-1].split('.')[0]
        # Initiate covariance model construction and calibration.
        cmc = CmConstructor(sto_path, model_out, name, cpu)
        # Construct the model.
        cmc.construct()
        # Calibrate the model.
        cmc.calibrate()
    print('\n### Construction of core set finished')

if __name__ == '__main__':
    main()
