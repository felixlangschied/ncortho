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
import pickle
import pyfaidx
import subprocess as sp
import sys

###############################################################################

# Parse a GTF file to store the coordinates for each protein-coding gene in a
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


# Perform reciprocal BLAST search and construct Stockholm alignment
def blastsearch(m_path, r_path, o_path, c, dust):

    core_set = {}

    with open(m_path) as infile:
        mirnas = [
            line.strip().split()
            for line in infile
            if not line.startswith('#')
        ]

    for mirna in mirnas:
        mirid = mirna[0]
        # Ensure that output folder exists and change to this folder.
        out_folder = '{}/{}'.format(o_path, mirid)
        if not os.path.isdir(out_folder):
            mkdir_cmd = 'mkdir -p {}'.format(out_folder)
            sp.call(mkdir_cmd, shell=True)
        os.chdir('{}/{}'.format(o_path, mirid))
        print('### {} ###'.format(mirid))
        # Coordinates of the ncRNA
        mchr = mirna[1]
        mstart = int(mirna[2])
        mend = int(mirna[3])

        # Start of reference bit score computation.
        # Convert RNA sequence into DNA sequence.
        preseq = mirna[5].replace('U', 'T')
        # print(preseq)
        # The miRNA precursors can show a low level of complexity, hence it is
        # required to deactivate the dust filter for the BLAST search.

        # check if blastdb of reference genome exists
        file_extensions = ['.nhr', '.nin', '.nsq']
        # file_extensions = ['.nhr']
        for fe in file_extensions:
            checkpath = '{}{}'.format(r_path, fe)
            if not os.path.isfile(checkpath):
                print(
                    'BLAST database does not exist for reference genome.\n'
                    'Constructing BLAST database.'
                )
                # At least one of the BLAST db files is not existent and has to be
                # created.
                db_command = (
                    'makeblastdb -dbtype nucl -in {}'
                        .format(r_path)
                )
                sp.call(db_command, shell=True)
                break
            else:
                print('BLAST database already exists.')

        bit_check = (
            'blastn -num_threads {0} -dust {1} -task megablast -db {2} '
            '-outfmt \"6 bitscore\"'.format(c, dust, r_path)
        )
        #print(bit_check)
        ref_bit_cmd = sp.Popen(
            bit_check, shell=True, stdin=sp.PIPE,
            stdout=sp.PIPE, stderr=sp.STDOUT, encoding='utf8'
        )
        ref_results, err = ref_bit_cmd.communicate(preseq)
        if not ref_results:
            print('# Sequence for {} not found in reference. Skipping..')
            skip_file = '{}/not_found_in_ref.fa'.format(o_path)
            with open(skip_file, 'w') as of:
                of.write('>{}\n{}\n'.format(mirid, preseq))
            continue
        ref_bit_score = float(ref_results.split('\n')[0].split('\t')[0])
        print(ref_bit_score)
        # End of reference bit score computation.
        print('Performing reciprocal BLAST search.')
        # 
        fasta = '{0}/{1}/{1}.fa'.format(o_path, mirid)
        #print(fasta)

        if os.path.isfile(fasta):
            print('FASTA file found for {}.'.format(mirna[0]))
            print('Checking BLAST database.')
# Check if BLAST database already exists, otherwise create it.
    # Database files are ".nhr", ".nin", ".nsq".
            file_extensions = ['.nhr', '.nin', '.nsq']
    #file_extensions = ['.nhr']
            for fe in file_extensions:
                checkpath = '{}{}'.format(fasta, fe)
                if not os.path.isfile(checkpath):
                    print(
                        'BLAST database does not exist.\n'
                        'Constructing BLAST database.'
                    )
            # At least one of the BLAST db files is not existent and has to be
            # created.
                    db_command = (
                        'makeblastdb -dbtype nucl -in {}'
                        .format(fasta)
                    )
                    sp.call(db_command, shell=True)
                    break
            else:
                print('BLAST database already exists.')

            blastn_cmd = (
                'blastn -num_threads {0} -task blastn -dust {1} -db {2} -outfmt \"6 '
                'sseqid evalue bitscore sseq\"'.format(c, dust, fasta)
            )
            blastn = sp.Popen(
                blastn_cmd, shell=True, stdin=sp.PIPE,
                stdout=sp.PIPE, stderr=sp.STDOUT, encoding='utf8'
            )
            results, err = blastn.communicate(preseq)
            # print(results)
##### Collect best hit for each core set species if it is within the accepted bit score range
            core_dict = {}
            result_list = results.split('\n')
            # print(result_list)
            if result_list:
                #print(result_list)
                for hit in result_list:
                    if hit:
                        hit_data = hit.split()
                        print(f'hit bitscore: {hit_data[2]}')
                        # print(f'you needed: {0.5*ref_bit_score}')
                        # print(f'and you got: {hit_data[2]}')
                        if (
                            not hit_data[0] in core_dict
                            and float(hit_data[2]) >= 0.5*ref_bit_score
                        ):
                        #if not hit_data[0] in core_dict:
                            core_dict[hit_data[0]] = hit
            if len(core_dict) == 0:
                print('# No region in the core species scored above the reference bit score threshold')
                #print(core_dict)
            #print(results.split('\n'))
            #print(results)
    ##### Re-BLAST #####
            #if core_dict:
            accept_dict = {}
            for species in core_dict.keys():
                print(species)
                # Make sure to eliminate gaps
                candidate_seq = core_dict[species].split()[3].replace('-', '')
                reblastn_cmd = (
                    'blastn -num_threads {0} -task blastn -dust {1} -db {2} -outfmt \"6'
                    ' sseqid sstart send evalue bitscore\"'
                    .format(c, dust, r_path)
                )
                reblastn = sp.Popen(
                    reblastn_cmd, shell=True, stdin=sp.PIPE,
                    stdout=sp.PIPE, stderr=sp.STDOUT, encoding='utf8'
                )
                reresults, reerr = reblastn.communicate(candidate_seq)
                print('Reverse search.')
                #print(reresults.split('\n')[0])
##### Check if reverse hit overlaps with reference miRNA
                if reresults:
                    first_hit = reresults.split('\n')[0].split()
                    rchr = first_hit[0]
                    rstart = int(first_hit[1])
                    rend = int(first_hit[2])
                    #print(first_hit[1], first_hit[2])
                    #print(first_hit)
                    #print('#####')
                    #print(rchr, mchr)
                    if rchr == mchr:
                        print('Same chromosome.')
                        if (
                            (rstart <= mstart and mstart <= rend)
                            or (rstart <= mend and mend <= rend)
                        ):
                    #first within second
                    #print("yes")
                    #print(tophit)
                            print('Reciprocity fulfilled.')
                            accept_dict[species] = core_dict[species]
                        elif (
                            (mstart <= rstart and rstart <= mend)
                            or (mstart <= rend and rend <= mend)
                        ):
                    #second within first
                    #print("yeah")
                    #print(tophit)
                            print('Reciprocity fulfilled.')
                            accept_dict[species] = core_dict[species]
                else:
                    # del core_dict[species]
                    print(
                        'No reverse hit for {}. Reciprocity unfulfilled.'
                        .format(mirid)
                    )

                #print(results.split('\n')[0].split())
            #for result in results.split('\n'):
            #    print(result)
        # If the FASTA file is not existent, the BLAST search cannot be
        # performed.
        else:
            print('FASTA file not found for {}.'.format(mirna[0]))
            continue
        corefile = '{}/{}_core.fa'.format(out_folder, mirid) 
        with open(corefile, 'w') as outfile:
            outfile.write('>{}\n{}\n'.format(mirid, preseq))
            print('##### Core set for {}: #####'.format(mirid))
            #print('>{}\n{}'.format(mirid, preseq))
            for accepted in accept_dict:
                #print('>{}\n{}'.format(accepted, core_dict[accepted].split('\t')[3]))
                outfile.write(
                    '>{}\n{}\n'
                    .format(accepted, accept_dict[accepted].split('\t')[3])
                    .replace('-', '')
                )
        alignment = '{}/{}.aln'.format(out_folder, mirid)
        stockholm = '{}/{}.sto'.format(out_folder, mirid)
        t_coffee = 't_coffee'
        #t_coffee = '/home/andreas/Applications/tcoffee/Version_11.00.8cbe486/bin/t_coffee'
        #t_coffee = '/home/andreas/Applications/T-COFFEE_installer_Version_12.00.7fb08c2_linux_x64/bin/t_coffee'
        #tc_cmd_1 = '{} -quiet=t_coffee.out -cpu 4 -special_mode=rcoffee -in {} -output=clustalw_aln > {}'.format(t_coffee, corefile, alignment)
        #tc_cmd_1 = '{} -quiet=t_coffee.out -cpu 4 -special_mode=rcoffee -in {} -output=clustalw_aln -outfile {}'.format(t_coffee, corefile, alignment)
        #tc_cmd_1 = '{} -quiet -multi_core={} -special_mode=rcoffee -in {} -output=clustalw_aln -outfile={}'.format(t_coffee, c, corefile, alignment)
        #tc_cmd_2 = '{} -other_pg seq_reformat -in {} -action +add_alifold -output stockholm_aln -out {}'.format(t_coffee, alignment, stockholm)
        # Call T-Coffee for the sequence alignment.
        print('Building T-Coffee alignment.')
        tc_cmd_1 = (
            '{} -quiet -multi_core={} -special_mode=rcoffee -in {} '
            '-output=clustalw_aln -outfile={}'
            .format(t_coffee, c, corefile, alignment)
        )
        sp.call(tc_cmd_1, shell=True)
        #with open(alignment) as aln_file:
            #aln_lines = [line for line in aln_file if not line.startswith('!')]
        #with open(alignment, 'w') as aln_file:
         #   for line in aln_lines:
          #      aln_file.write(line)
        # Extend the sequence-based alignment by structural information.
        # Create Stockholm alignment.
        print('Adding secondary structure to Stockholm format.')
        tc_cmd_2 = (
            '{} -other_pg seq_reformat -in {} -action +add_alifold -output '
            'stockholm_aln -out {}'
            .format(t_coffee, alignment, stockholm)
        )
        sp.call(tc_cmd_2, shell=True)


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
    # Maximum gene insertions
    parser.add_argument(
        '--dust', metavar='yes/no', type=str,
        help='Use BLASTn dust filter. Decreases number of models created but improves runtime and possibly specificity', nargs='?',
        const='yes', default='yes'
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
    ref_gtf_path = args.reference
    ref_genome = args.genome
    oma_paths = glob.glob('{}/*'.format(args.pairwise))
    core_gtf_paths = args.core
    core_fa_paths = args.query
    mgi = args.mgi
    dust = args.dust

###############################################################################

# Parse the pairwise orthologs
    for oma_path in oma_paths:
        taxon = oma_path.split('/')[-1]
        with open(oma_path) as omafile:
#        orthologs = {ref: core for (ref, core) in 
#[(line.split()[0], line.split()[1]) for line in omafile.readlines()
#if len(line.split()) == 4]}
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
    ref_dict = gtf_parser(ref_gtf_path)

#Determine the position of each miRNA and its neighboring gene(s)
    for mirna in mirnas:
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
                        print(left_hits)
                        print(right_hits)
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

    print(neighbor_dict)

### Search for the coordinates of the orthologs and extract the sequences
    print('# starting now with coordinate search\n')
    for taxon in neighbor_dict:
        print('Starting synteny analysis for {}'.format(taxon))
        gtf_path = '{0}/{1}.gtf'.format(core_gtf_paths, taxon)
        fasta_path = glob.glob('{0}/{1}*.fa'.format(core_fa_paths, taxon))
        #print(gtf_path)
        #print(fasta_path)
        if len(fasta_path) != 1:
            print('Unable to identify genome file for {}'.format(taxon))
            sys.exit()
        genome = pyfaidx.Fasta(fasta_path[0])
        print('Trying to parse GTF file for {}.'.format(taxon))
        try:
            core_gtf_dict = gtf_parser(gtf_path)
            print('Done')
        except:
            print('No GTF file found for {}'.format(taxon))
            sys.exit()
        # print('YOU MADE IT THIS FAR.')
        for mirna in neighbor_dict[taxon]:
            print(mirna)
            style = neighbor_dict[taxon][mirna][0]
            print(f'# style = {style}')
            if style == 'inside' or style == 'opposite':
###############################################################################
                try:
                    ortho_data = (
                        core_gtf_dict[neighbor_dict[taxon][mirna][1]]
                    )
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
                    except:
                        mirna_dict[mirna] = {taxon: seq}
                except:
                    print('{} not found in GTF file.'.format(mirna[1]))
            elif style == 'in-between':
                left_data = (
                    core_gtf_dict[neighbor_dict[taxon][mirna][1][0]]
                )
                right_data = (
                    core_gtf_dict[neighbor_dict[taxon][mirna][1][1]]
                )
                print('#########################')
# Test to see if the two orthologs are themselves neighbors where their
# distance cannot be larger than the selected mgi value. This accounts
# for insertions in the core species.
# TODO: Apply mgi also to the reference species to account for insertions
# in the reference.
                print(f'left_data: {left_data}')
                print(f'right_data: {right_data}')
                if (
                    left_data[0] == right_data[0]
                    and abs(left_data[1] - right_data[1]) <= mgi
                ):
# Determine which sequence to include for the synteny-based ortholog search
# depending on the order of orthologs. The order of the orthologs in the core
# species might be inverted compared to that in the reference species.
###############################################################################
                    if left_data[1] < right_data[1]:
                        print('left')
                        print(core_gtf_dict[left_data[0]][left_data[1]])
                        print(core_gtf_dict[right_data[0]][right_data[1]])
                        contig = left_data[0]
                        print(contig)
                        seq_start = int(
                            core_gtf_dict[left_data[0]][left_data[1]][2]
                        )
                        print(seq_start)
                        seq_end = (
                            core_gtf_dict[right_data[0]][right_data[1]][1]
                        )
                        print(seq_end)
                        seq = genome[contig][seq_start-1:seq_end].seq
                        try:
                            mirna_dict[mirna][taxon] = seq
                        except:
                            mirna_dict[mirna] = {taxon: seq}
                    elif right_data[1] < left_data[1]:
                        print('right')
                        print(core_gtf_dict[left_data[0]][left_data[1]])
                        print(core_gtf_dict[right_data[0]][right_data[1]])
                        contig = left_data[0]
                        print(contig)
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
                        print(seq_start)
                        seq_end = (
                            core_gtf_dict[left_data[0]][left_data[1]][1]
                        )
                        print(seq_end)
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
                    print(left_data)
                    print(right_data)
            else:
                print('## Neither inside, opposite, nor in-between')
                print(neighbor_dict[taxon][mirna])
                print('##')

            # continue
    # print(mirna_dict.keys())
    #print(mirna_dict['mmu-mir-1224'])
    # print(mirna_dict['mmu-mir-6899'].keys())

    def write_output():
        for mirna in mirna_dict:
            with open('{0}/{1}/{1}.fa'.format(output, mirna), 'w') as outfile:
                for core_taxon in mirna_dict[mirna]:
                    outfile.write(
                        '>{0}\n{1}\n'
                        .format(core_taxon, mirna_dict[mirna][core_taxon])
                    )

    write_output()

    #fasta_path = '/home/andreas/Documents/Internship/ncOrtho_to_distribute/ncortho_python/test_core_set_construction/test_fasta'
    #blastsearch(mirna_path, fasta_path, ref_genome, output, cpu)
    print('# Starting reciprocal blast\n')
    blastsearch(mirna_path, ref_genome, output, cpu, dust)

if __name__ == '__main__':
    main()
