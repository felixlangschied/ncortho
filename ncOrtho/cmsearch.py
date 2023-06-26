import pyfaidx
import subprocess as sp
import os
import logging


def heuristic_search(blastdb, rna, query, heuristic, cpu, out):

    def candidate_blast(db, c, evalue, seq, task):
        # extract candidate regions
        blast_command = (
            'blastn -task {0} -db {1} '
            '-num_threads {2} -evalue {3} '
            '-outfmt "6 sseqid sstart send sstrand length"'.format(task, db, c, evalue)
        )
        blast_call = sp.Popen(
            blast_command, shell=True, stdin=sp.PIPE, stdout=sp.PIPE, stderr=sp.PIPE, encoding='utf8'
        )
        res, err = blast_call.communicate(seq)
        if err:
            raise sp.SubprocessError(err)

        result = list(filter(None, res.split('\n')))
        hit_list = [line.split() for line in result if line and float(line.split()[-1]) >= blast_len_cut]

        return hit_list

    def extract_candidate_regions(candidate_list, query, extraregion=1000):
        hit_at_start = False
        hit_at_end = False
        genome = pyfaidx.Fasta(query)

        regions = []
        for hit in candidate_list:
            chrom, start, end, strand, length = hit
            strand = strand.replace('plus', '+').replace('minus', '-')
            if strand == '-':
                start, end = end, start
            start = int(start)
            end = int(end)

            if start > extraregion:
                n_start = start - extraregion
            else:
                # hit starts at the beginning of the chromosome
                n_start = 0
                hit_at_start = True
            n_end = end + extraregion
            if n_end > genome[chrom][-1].end:
                # hit is at the end of the chromosome
                n_end = genome[chrom][-1].end
                hit_at_end = True
            if strand == '+':
                sequence = genome[chrom][n_start:n_end].seq
            elif strand == '-':
                sequence = genome[chrom][n_start:n_end].reverse.complement.seq
            else:
                raise ValueError(f'Unknown strand {strand}')
            regions.append(f">{chrom}|{str(start)}|{str(end)}|{strand}|{hit_at_start}|{hit_at_end}\n{sequence}\n")

        return regions


    ################################################################################################
    blast_len_cut = len(rna.seq) * heuristic[2]
    candidate_list = candidate_blast(blastdb, cpu, heuristic[1], rna.seq, 'blastn')
    if not candidate_list:
        return ''

    candidate_regions = extract_candidate_regions(candidate_list, query)
    fasta = f'{out}/candidate_regions_{rna.name}.fa'
    with open(fasta, 'w') as of:
        for line in candidate_regions:
            of.write(line)

    return fasta


#[
# ['(1)', '!', '3.7e-23', '94.0', '0.0', 'NW_020173731.1|9158415|9158472|+|False|False', '1000', '1057', '+', 'cm', 'no', '0.41', '-'],
# ['(2)', '!', '3.7e-23', '94.0', '0.0', 'NW_020173731.1|9274958|9275015|+|False|False', '1000', '1057', '+', 'cm', 'no', '0.41', '-'],
# ['(3)', '!', '8.3e-15', '62.1', '0.0', 'NW_020172457.1|4920040|4920099|+|False|False', '1000', '1059', '+', 'cm', 'no', '0.33', '-']
# ]

def perform_model_search(models, rna, candidate_region_fasta, output, cm_cutoff, cpu, phmm):

    def parse_phmm(tblout):
        phmm_results = []
        with open(tblout) as fh:
            for line in fh:
                if line.startswith('#'):
                    continue
                data = line.strip().split()
                # print(data)
                info = data[0]
                start, end = data[6:8]
                strand, evalue, score = data[11:14]
                phmm_results.append(['a', 'b', evalue, score, 'c', info, start, end, strand, 'd', 'e', 'f', 'g'])
        return phmm_results

    def parse_cm(res):
        cm_results = []
        for line in res.split('\n'):
            # print(line)
            if line.startswith('  ('):
                if 'No hits detected that satisfy reporting thresholds' in line:
                    return []
                else:
                    data = line.strip().split()
                    cm_results.append(data)
        return cm_results

    ################################################################################
    cut_off = rna.bit * cm_cutoff

    if phmm:
        tblout = f'{output}/{rna.name}_phmm.out'
        call = sp.Popen(
            f'nhmmer --cpu {cpu} --tblout {tblout} -T {cut_off} --incT {cut_off} {models}/{rna.name}.phmm {candidate_region_fasta}',
            shell=True, stdin=sp.PIPE, stdout=sp.PIPE, stderr=sp.PIPE, encoding='utf8'
        )
        res, err = call.communicate(rna.seq)
        if err:
            raise sp.SubprocessError(err)
        results = parse_phmm(tblout)
        os.remove(tblout)
    else:
        call = sp.Popen(
            f'cmsearch --cpu {cpu} --noali -T {cut_off} --incT {cut_off} {models}/{rna.name}.cm {candidate_region_fasta}',
            shell=True, stdout=sp.PIPE, stderr=sp.PIPE, encoding='utf8'
        )
        res, err = call.communicate(rna.seq)
        if err:
            raise sp.SubprocessError(err)
        results = parse_cm(res)

    return results


def parse_cmsearch(cm_results):
    logger = logging.getLogger('ncortho')
    logger.setLevel(level=logging.DEBUG)
    return_data = []
    for data in cm_results:
        score = float(data[3])
        chrom = data[5]
        start, end = map(int, data[6:8])
        strand = data[8]
        return_data.append([chrom, start, end, strand, score])
    return return_data


def parse_cmsearch_for_heuristic(cm_results, extraregion=1000):
    logger = logging.getLogger('ncortho')
    logger.setLevel(level=logging.DEBUG)
    return_data = []
    for data in cm_results:
        # parse CMsearch output of candidate regions to acutal genomic regions
        blast_chrom = data[5].split('|')[0]
        blast_start, blast_end = map(int, data[5].split('|')[1:3])
        hit_at_start, hit_at_end = data[5].split('|')[4:6]
        blast_strand = data[5].split('|')[3]
        cm_start, cm_end = map(int, data[6:8])
        # cm_start, cm_end = map(int, data[5:7])
        cm_strand = data[8]
        if cm_strand == '-':
            # reverse hits should not be possible since blasthits were extracted
            # as reverse complement if they were on the minus strand
            # logger.info('Unexpected hit of CMsearch')
            continue
        if hit_at_start == 'True':
            # print(f'# cmsearch hit for {rna_id} at the beginning of a chromosome in the query species')
            hit_start = cm_start
            hit_end = cm_end
        elif hit_at_end == 'True':
            # print(f'# cmsearch hit for {rna_id} at the end of a chromosome in the query species')
            hit_start = blast_start + (cm_start - extraregion) - 1
            hit_end = hit_start + (cm_end - cm_start) - 1
        else:
            hit_start = blast_start + (cm_start - extraregion) - 1
            hit_end = hit_start + (cm_end - extraregion) - 1
        # fill output variable
        return_data.append([blast_chrom, hit_start, hit_end, blast_strand, float(data[3])])
    return return_data


def remove_duplicates_apply_cutoffs(cms, cmc, lc, rna):
    """
    Parse the output of cmsearch while eliminating duplicates and filtering
    entries according to the defined cutoff.

    Parameters
    ----------
    cms     :   List with CMsearch hits [chrom, start, end, strand, score]
    cmc     :   Cutoff to decide which candidate hits should be included for the reverse BLAST search
    lc      :   Length cutoff
    rna     :   rna class

    Returns
    -------
    hits_dict : Dictionary containing hits of the cmsearch

    """
    # Output
    hits_dict = {}
    # Required for finding duplicates, stores hits per chromosome
    chromo_dict = {}

    for candidate_nr, hit in enumerate(cms, 1):
        # print(hit)
        h_chrom, h_start, h_end, h_strand, h_score = hit
        # blastparser expects the start to be the smaller number, will extract reverse complement if on - strand
        if h_start > h_end:
            h_start, h_end = h_end, h_start
        length = h_end - h_start

        if length <= length * lc:
            continue
        if h_score <= rna.bit * cmc:
            continue

        candidate = f'{rna.name}_c{candidate_nr}'
        data = candidate, h_chrom, h_start, h_end, h_strand, h_score
        hits_dict[candidate] = data

        # Store the hits that satisfy the bit score cutoff to filter
        # duplicates.
        if h_chrom in chromo_dict:
            chromo_dict[h_chrom] = []
        chromo_dict[h_chrom] = [data]
    # print(hits_dict)

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
                                del hits_dict[chromo_dict[chromo][chitnr][0]]

                            else:
                                del hits_dict[chromo_dict[chromo][hitnr][0]]

    return hits_dict


def model_search(rna, cm_cutoff, cpu, msl, models, query, blastdb, out, cleanup, heuristic_col, phmm):
    """
    Parameters
    ----------
    rna       :   rna object
    cm_cutoff   :   Minimum bit score of the CMsearch hit relative to the
                    bit score that results from searching in the reference genome with the same model
    cpu         :   Number of cores to use
    msl         :   Minimum length of CMsearch hit relative to the reference pre-rna length
    models      :   Path to the directory that contains the CMs
    query       :   Path to the query genome
    blastdb     :   Optional Path to a BLASTdb of the query genome
    out         :   Path to the output directory
    cleanup     :   Delete intermediate files True/False
    heuristic   :   Should heuristic mode be used and if yes, set parameters (True/False, evalue, minlength)

    Returns
    -------
    cm_results : Dictionary containing hits of the cmsearch

    """
    # print(os.path.isfile(f'{models}/{rna.name}.phmm'))
    if (
            not phmm and not os.path.isfile(f'{models}/{rna.name}.cm') or
            phmm and not os.path.isfile(f'{models}/{rna.name}.phmm')
    ):
        return False, 'No model found'
    heuristic, heuristic_evalue, heuristic_length, sensitive_heuristic = heuristic_col
    if heuristic:
        candidate_region_fasta = heuristic_search(blastdb, rna, query, heuristic_col, cpu, out)

        if not candidate_region_fasta and sensitive_heuristic:
            # No candidate regions above length cutoff found with BLASTn. Trying again without heuristic.
            cm_results = perform_model_search(models, rna, query, out, cm_cutoff, cpu, phmm)
        elif not candidate_region_fasta:
            return False, 'No candidate region found with BLASTn'
        else:
            cm_results = perform_model_search(models, rna, candidate_region_fasta, out, cm_cutoff, cpu, phmm)
            cm_results = parse_cmsearch_for_heuristic(cm_results)

        if cleanup:
            os.remove(candidate_region_fasta)
    else:
        cm_results = perform_model_search(models, rna, query, out, cm_cutoff, cpu, phmm)
        cm_results = parse_cmsearch(cm_results)

    if not cm_results:
        return False, 'No CMsearch results in candidate regions or query genome'

    parsed_cm_results = remove_duplicates_apply_cutoffs(cm_results, cm_cutoff, msl, rna)
    if not parsed_cm_results:
        return False, 'No CMsearch results above threshold'
    else:
        return parsed_cm_results, 'Sucess'










