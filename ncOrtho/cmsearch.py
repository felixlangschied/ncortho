import glob
import pyfaidx
import subprocess as sp
import sys
import os


def cmsearcher(mirna, cm_cutoff, cpu, msl, models, query, out, cleanup, heuristic):
    mirna_id = mirna.name

    blastdb = '{}/query_BLASTdb/{}'.format(out, query.split('/')[-1])
    # Calculate the bit score cutoff.
    cut_off = mirna.bit * cm_cutoff
    # Calculate the length cutoff.
    len_cut = len(mirna.pre) * msl
    # Perform covariance model search.
    # Report and inclusion thresholds set according to cutoff.

    if heuristic:
        # Test if BLASTdb exists for query species
        file_extensions = ['nhr', 'nin', 'nsq']
        for fe in file_extensions:
            files = glob.glob(f'{blastdb}*{fe}')
            if not files:
                # At least one of the BLAST db files is not existent and has to be
                # created.
                db_command = 'makeblastdb -in {} -out {} -dbtype nucl'.format(query, blastdb)
                sp.call(db_command, shell=True)
                break
        # extract candidate regions
        print('# Identifying candidate regions for cmsearch heuristic')
        tmp_out = '{0}/tmp_blast_{1}.out'.format(out, mirna_id)
        with open(tmp_out, 'w') as inf:
            inf.write(">{}\n{}\n".format(mirna_id, mirna.pre))
        blast_command = (
            'blastn -task blastn -db {0} -query {1} '
            '-num_threads {2} -outfmt "6 sseqid sstart send sstrand"'.format(blastdb, tmp_out, cpu)
        )
        res = sp.run(blast_command, shell=True, capture_output=True)
        os.remove(tmp_out)
        if res.stdout:
            print('Blast step finished')
        else:
            cm_results = False
            return cm_results
        result = res.stdout.decode('utf-8').split('\n')
        hit_list = [line.split() for line in result if line]
        # open indexed query fasta
        genes = pyfaidx.Fasta(query)
        # collect heuristic sequences
        print('collecting sequences')
        heuristic_fa = '{0}/heuristic_{1}.out'.format(out, mirna_id)
        with open(heuristic_fa, 'w') as of:
            for hit in hit_list:
                chrom, start, end, strand = hit
                strand = strand.replace('plus', '+')
                strand = strand.replace('minus', '-')
                header = ">{}|{}|{}|{}\n".format(chrom, start, end, strand)
                # change start and end point to include neighbourhood
                if strand == '+':
                    n_start = int(start) - 1000
                    n_end = int(end) + 1000
                    sequence = genes[chrom][n_start:n_end].seq
                elif strand == '-':
                    n_start = int(start) + 1000 
                    n_end = int(end) - 1000
                    sequence = genes[chrom][n_end:n_start].seq#.reverse.complement.seq
                else:
                    print('Something went wrong during the heuristic cmsearch')
                    sys.exit()
                of.write(header)
                of.write(f'{sequence}\n')
        # start the cmsearch with the heuristic candidates
        heuristic_cms = '{0}/cmsearch_heuristic.out'.format(out, mirna_id)
        # if not os.path.isfile(cms_output):
        print('# Running covariance model search for {}'.format(mirna_id))
        cms_command = (
            'cmsearch -T {5} --incT {5} --cpu {0} --noali '
            '--tblout {1} {2}/{3}.cm {4}'
                .format(cpu, heuristic_cms, models, mirna_id, heuristic_fa, cut_off)
        )
        res = sp.run(cms_command, shell=True, capture_output=True)
        if not res.stdout:
            print('# cmsearch did not find anything')
            cm_results = False
            return cm_results

        # else:
        #     print('# Found cm_search results at: {}. Using those'.format(cms_output))
        #     cm_results = cmsearch_parser(cms_output, cut_off, len_cut, mirna_id)
        heur_results = heuristic_cms.replace('.out', '_res.out')
        with open(heuristic_cms, 'r') as inf, open(heur_results, 'w') as of:
            for line in inf:
                if not line.startswith('#'):
                    data = line.strip().split()
                    blast_chrom = data[0].split('|')[0]
                    blast_start, blast_end = map(int, data[0].split('|')[1:3])
                    # blast_strand = data[0].split('|')[-1]
                    cm_start, cm_end = map(int, data[7:9])
                    # cm_strand = data[9]
                    # if not cm_strand == blast_strand:
                    #     print('Rejecting cmsearch hit. Not on same strand as BLAST hit')
                    #     break
                    hit_start = blast_start + (cm_start - 1000)
                    hit_end = blast_end + (cm_end - 1000)
                    data[0] = blast_chrom
                    data[7] = hit_start
                    data[8] = hit_end
                    data = [str(entry) for entry in data]
                    of.write('\t'.join(data))
                    of.write('\n')
        cm_results = cmsearch_parser(heur_results, cut_off, len_cut, mirna_id)
        if cleanup:
            os.remove(heuristic_fa)
            os.remove(heur_results)
            os.remove(heuristic_cms)
    else:
        cms_output = '{0}/cmsearch_{1}.out'.format(out, mirna_id)
        if not os.path.isfile(cms_output):
            print('# # Running covariance model search for {}'.format(mirna_id))
            cms_command = (
                'cmsearch -T {5} --incT {5} --cpu {0} --noali '
                '--tblout {1} {2}/{3}.cm {4}'
                    .format(cpu, cms_output, models, mirna_id, query, cut_off)
            )
            sp.call(cms_command, shell=True)
        else:
            print('# Found cm_search results at: {}. Using those'.format(cms_output))
        cm_results = cmsearch_parser(cms_output, cut_off, len_cut, mirna_id)

    return cm_results

# cmsearch_parser: Parse the output of cmsearch while eliminating
#                  duplicates and filtering entries according to the
#                  defined cutoff.
# Arguments:
# cms: path to cmsearch output
# cmc: cutoff to decide which candidate hits should be included for the
#      reverse BLAST search
# lc: length cutoff
# mirid: name/id of the microRNA
def cmsearch_parser(cms, cmc, lc, mirid):
    """
    Parse the output of cmsearch while eliminating duplicates and filtering
    entries according to the defined cutoff.

    Parameters
    ----------
    cms : TYPE
        DESCRIPTION.
    cmc : TYPE
        DESCRIPTION.
    lc : TYPE
        DESCRIPTION.
    mirid : TYPE
        DESCRIPTION.

    Returns
    -------
    hits_dict : TYPE
        DESCRIPTION.

    """
    # Output
    hits_dict = {}
    # Required for finding duplicates, stores hits per chromosome
    chromo_dict = {}
    cut_off = cmc

    with open(cms) as cmsfile:
        # Collect only the hits which satisfy the bit score cutoff.
        hits = [
            line.strip().split() for line in cmsfile
            if not line.startswith('#')
               and float(line.strip().split()[14]) >= cut_off
               and abs(int(line.split()[7]) - int(line.strip().split()[8])) >= lc
        ]
        # Add the hits to the return dictionary.
        if hits:
            for candidate_nr, hit in enumerate(hits, 1):
                data = (
                    '{0}_c{1}'.format(mirid, candidate_nr),
                    hit[0], hit[7], hit[8], hit[9], hit[14]
                )
                hits_dict[data[0]] = data

                # Store the hits that satisfy the bit score cutoff to filter
                # duplicates.
                try:
                    chromo_dict[data[1]].append(data)
                except:
                    chromo_dict[data[1]] = [data]

    # Loop over the candidate hits to eliminate duplicates.
    for chromo in chromo_dict:
        nrhits = len(chromo_dict[chromo])
        if nrhits > 1:
            for hitnr in range(nrhits):
                start = int(chromo_dict[chromo][hitnr][2])
                stop = int(chromo_dict[chromo][hitnr][3])
                strand = chromo_dict[chromo][hitnr][4]
                score = float(chromo_dict[chromo][hitnr][5])
                for chitnr in range(hitnr +1, nrhits):
                    if strand != chromo_dict[chromo][chitnr][4]:
                        cstart = int(chromo_dict[chromo][chitnr][2])
                        cstop = int(chromo_dict[chromo][chitnr][3])
                        cscore = float(chromo_dict[chromo][chitnr][5])
                        # Test if the two hits from opposite strands
                        # overlap, which means one of them is
                        # (probably) a false positive.
                        # Out of two conflicting hits, the one with the
                        # highest cmsearch bit score is retained.
                        if (
                                start in range(cstart, cstop +1)
                                or stop in range(cstart, cstop +1)
                                or cstart in range(start, stop +1)
                                or cstop in range(start, stop +1)
                        ):
                            if score > cscore:
                                try:
                                    del hits_dict[chromo_dict[chromo]
                                    [chitnr][0]]
                                except:
                                    pass
                            else:
                                try:
                                    del hits_dict[chromo_dict[chromo]
                                    [hitnr][0]]
                                except:
                                    pass
    return hits_dict
