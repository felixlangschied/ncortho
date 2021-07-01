import subprocess as sp
import os
import glob
import sys


# Perform reciprocal BLAST search and construct Stockholm alignment
def blastsearch(mirna, r_path, o_path, c, dust):
    """

    Parameters
    ----------
    m_path  :   Path to file with miRNA information
    r_path  :   Path to reference genome.
    o_path  :   Output Name
    c       :   Number of Threads for BLAST search.
    dust    :   Use dust filter for BLAST searches (yes/no)

    Returns
    -------
    Path to Alignment of core pre-miRNAs in Stockholm format

    """
    core_set = {}


    # with open(m_path) as infile:
    #     mirnas = [
    #         line.strip().split()
    #         for line in infile
    #         if not line.startswith('#')
    #     ]
    # for mirna in mirnas:
    mirid = mirna[0]
    # this fasta file is created by the the main() script
    fasta = '{0}/{1}/{1}.fa'.format(o_path, mirid)
    if not os.path.isfile(fasta):
        print('ERROR: FASTA file not found for {}.'.format(mirna[0]))
        return None

    # Ensure that output folder exists and change to this folder.
    out_folder = '{}/{}'.format(o_path, mirid)
    if not os.path.isdir(out_folder):
        mkdir_cmd = 'mkdir -p {}'.format(out_folder)
        sp.call(mkdir_cmd, shell=True)
    os.chdir('{}/{}'.format(o_path, mirid))
    print('\n### {} ###'.format(mirid))
    # Coordinates of the ncRNA
    mchr = mirna[1].replace('chr', '')
    mstart = int(mirna[2])
    mend = int(mirna[3])

    # Start of reference bit score computation.
    # Convert RNA sequence into DNA sequence.
    preseq = mirna[5].replace('U', 'T')

    # The miRNA precursors can show a low level of complexity, hence it is
    # required to deactivate the dust filter for the BLAST search.


    # check if blastdb of reference genome exists
    fname = '.'.join(r_path.split("/")[-1].split('.')[0:-1])
    ref_blastdb = f'{o_path}/refBLASTdb/{fname}'
    file_extensions = ['.nhr', '.nin', '.nsq']
    count = 0
    for fe in file_extensions:
        files = glob.glob(f'{ref_blastdb}*{fe}')
        if not files:
            print(
                'BLAST database does not exist for reference genome.\n'
                'Constructing BLAST database.'
            )
            # At least one of the BLAST db files is not existent and has to be
            # created.
            db_command = 'makeblastdb -in {} -out {} -dbtype nucl'.format(r_path, ref_blastdb)
            sp.call(db_command, shell=True)
            break
        else:
            count += 1
    if count == 3:
        print('BLAST database of reference genome already exists.')

    bit_check = (
        'blastn -num_threads {0} -dust {1} -task megablast -db {2} '
        '-outfmt \"6 bitscore\"'.format(c, dust, ref_blastdb)
    )
    print('# Determining reference bit score..')
    ref_bit_cmd = sp.Popen(
        bit_check, shell=True, stdin=sp.PIPE,
        stdout=sp.PIPE, stderr=sp.STDOUT, encoding='utf8'
    )
    ref_results, err = ref_bit_cmd.communicate(preseq)
    if err:
        print(err)
    if not ref_results:
        print('WARNING: Reference sequence of {} not found in reference Genome. '
              'Consider turning the dust filter off'.format(mirid))
        skip_file = '{}/not_found_in_ref.fa'.format(o_path)
        if not os.path.isfile(skip_file):
            with open(skip_file, 'w') as of:
                of.write('>{}\n{}\n'.format(mirid, preseq))
        else:
            with open(skip_file, 'a') as of:
                of.write('>{}\n{}\n'.format(mirid, preseq))
        return None
    ref_bit_score = float(ref_results.split('\n')[0].split('\t')[0])
    print(f'BLAST bit score of reference miRNA against reference genome: {ref_bit_score}')

    print('# Trying to find pre-miRNA candidates in syntenic region')
    # Check if BLAST database already exists, otherwise create it.
    file_extensions = ['.nhr', '.nin', '.nsq']
    for fe in file_extensions:
        checkpath = '{}{}'.format(fasta, fe)
        if not os.path.isfile(checkpath):
            print(
                'Constructing BLAST database of candidate regions.'
            )
    # At least one of the BLAST db files is not existent and has to be
    # created.
            db_command = (
                'makeblastdb -dbtype nucl -in {}'
                .format(fasta)
            )
            sp.call(db_command, shell=True)
            break

    blastn_cmd = (
        'blastn -num_threads {0} -task blastn -dust {1} -db {2} -outfmt \"6 '
        'sseqid evalue bitscore sseq\"'.format(c, dust, fasta)
    )
    blastn = sp.Popen(
        blastn_cmd, shell=True, stdin=sp.PIPE,
        stdout=sp.PIPE, stderr=sp.PIPE, encoding='utf8'
    )
    results, err = blastn.communicate(preseq)
    if err:
        print(err)
        sys.exit()


    ##### Collect best hit for each core set species if it is within the accepted bit score range
    core_dict = {}
    result_list = results.split('\n')
    # print(result_list)
    if result_list:
        #print(result_list)
        for hit in result_list:
            if hit:
                hit_data = hit.split()
                # print(hit_data)
                # print(f'hit bitscore: {hit_data[2]}')
                # print(f'you needed: {0.5*ref_bit_score}')
                # print(f'and you got: {hit_data[2]}')
                if (
                    not hit_data[0] in core_dict
                    and float(hit_data[2]) >= 0.5*ref_bit_score
                ):
                    core_dict[hit_data[0]] = hit
                    print(f'pre-miRNA candidate found for {hit_data[0]}! ({hit_data[2]}/{0.5*ref_bit_score})')
                    break
                else:
                    print(f'Syntenic candidate region BLAST hit below reference bit score threshold ({hit_data[2]}/{0.5*ref_bit_score})')

    if len(core_dict) == 0:
        print('WARNING: No region in the core species scored above the reference bit score threshold')


    ##### Re-BLAST #####
    print('\n### Starting reciprocal BLAST search.')
    accept_dict = {}
    for species in core_dict.keys():
        print(f'# {species}')
        # Make sure to eliminate gaps
        candidate_seq = core_dict[species].split()[3].replace('-', '')
        reblastn_cmd = (
            'blastn -num_threads {0} -task blastn -dust {1} -db {2} -outfmt \"6'
            ' sseqid sstart send evalue bitscore\"'
            .format(c, dust, ref_blastdb)
        )
        reblastn = sp.Popen(
            reblastn_cmd, shell=True, stdin=sp.PIPE,
            stdout=sp.PIPE, stderr=sp.STDOUT, encoding='utf8'
        )
        reresults, reerr = reblastn.communicate(candidate_seq)
        if reerr:
            print(reerr)
            sys.exit()

    # Check if reverse hit overlaps with reference miRNA
        if reresults:
            first_hit = reresults.split('\n')[0].split()
            rchr = first_hit[0]
            rstart = int(first_hit[1])
            rend = int(first_hit[2])
            if rchr == mchr:
                # print('Same chromosome.')
                if (
                    (rstart <= mstart and mstart <= rend)
                    or (rstart <= mend and mend <= rend)
                ):
                    #first within second
                    print('Reciprocity fulfilled.')
                    accept_dict[species] = core_dict[species]
                elif (
                    (mstart <= rstart and rstart <= mend)
                    or (mstart <= rend and rend <= mend)
                ):

                    # second within first
                    print('Reciprocity fulfilled.')
                    accept_dict[species] = core_dict[species]
                else:
                    print('Reciprocity unfulfilled.')
            else:
                print('Hit on chromosome {}, expected {}'.format(rchr, mchr))
        else:
            print(
                'No reverse hit for {}. Reciprocity unfulfilled.'
                .format(mirid)
            )
    # write intermediate output file
    corefile = '{}/{}_core.fa'.format(out_folder, mirid)
    with open(corefile, 'w') as outfile:
        outfile.write('>{}\n{}\n'.format(mirid, preseq))
        print('\n### Starting Alignment')
        # print('>{}\n{}'.format(mirid, preseq))
        for accepted in accept_dict:
            # print('>{}\n{}'.format(accepted, core_dict[accepted].split('\t')[3]))
            outfile.write(
                '>{}\n{}\n'
                    .format(accepted, accept_dict[accepted].split('\t')[3])
                    .replace('-', '')
            )
    alignment = '{}/{}.aln'.format(out_folder, mirid)
    stockholm = '{}/{}.sto'.format(out_folder, mirid)
    t_coffee = 't_coffee'

    # Call T-Coffee for the sequence alignment.
    print('Building T-Coffee alignment.')
    tc_cmd_1 = (
        '{} -quiet -multi_core={} -special_mode=rcoffee -in {} '
        '-output=clustalw_aln -outfile={}'
            .format(t_coffee, c, corefile, alignment)
    )
    sp.call(tc_cmd_1, shell=True)

    # Extend the sequence-based alignment by structural information.
    # Create Stockholm alignment.
    print('Adding secondary structure to Stockholm format.')
    tc_cmd_2 = (
        '{} -other_pg seq_reformat -in {} -action +add_alifold -output '
        'stockholm_aln -out {}'
            .format(t_coffee, alignment, stockholm)
    )
    sp.call(tc_cmd_2, shell=True)
    return stockholm


