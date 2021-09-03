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

import pyfaidx
import subprocess as sp
import os
import sys


def candidate_blast(db, c, evalue, seq, task):
    # extract candidate regions
    print('# Identifying candidate regions for cmsearch heuristic')
    blast_command = (
        'blastn -task {0} -db {1} '
        '-num_threads {2} -evalue {3} '
        '-outfmt "6 sseqid sstart send sstrand length"'.format(task, db, c, evalue)
    )
    blast_call = sp.Popen(
        blast_command, shell=True, stdin=sp.PIPE, stdout=sp.PIPE, stderr=sp.PIPE, encoding='utf8'
    )
    res, err = blast_call.communicate(seq)
    return res, err



def cmsearcher(mirna, cm_cutoff, cpu, msl, models, query, blastdb, out, cleanup, heuristic):
    """
    Parameters
    ----------
    mirna       :   Mirna object
    cm_cutoff   :   Minimum bit score of the CMsearch hit relative to the
                    bit score that results from searching in the reference genome with the same model
    cpu         :   Number of cores to use
    msl         :   Minimum length of CMsearch hit relative to the reference pre-miRNA length
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
    mirna_id = mirna.name
    cm_results = False
    # Calculate the bit score cutoff.
    cut_off = mirna.bit * cm_cutoff
    # Calculate the length cutoff.
    len_cut = len(mirna.pre) * msl
    blast_len_cut = len(mirna.pre) * heuristic[2]
    # if no model exists for mirna return False
    if not os.path.isfile('{}/{}.cm'.format(models, mirna_id)):
        print('# No model found for {}. Skipping..'.format(mirna_id))
        exitstatus = 'No covariance model found'
        return cm_results, exitstatus
    # Perform covariance model search.
    # Report and inclusion thresholds set according to cutoff.

    if heuristic[0]:
        # extract candidate regions
        #candidate_blast(db, c, evalue, seq, task):
        res, err = candidate_blast(blastdb, cpu, heuristic[1], mirna.pre, 'blastn')
        if not err and not res:
            print('No BLASTn candidate regions with pre-miRNA. Trying mature sequence')
            res, err = candidate_blast(blastdb, cpu, 10, mirna.mature, 'blastn-short')
            # turn off length blast length cutoff
            blast_len_cut = 0
        if err:
            print(f'ERROR: {err}')
            exitstatus = 'ERROR during candidate region BLAST search'
            return cm_results, exitstatus

        # print('Blast step finished')
        result = list(filter(None, res.split('\n')))
        hit_list = [line.split() for line in result if line and float(line.split()[-1]) >= blast_len_cut]
        if not hit_list:
            exitstatus = 'No BLASTn candidate regions above length cutoff'
            return cm_results, exitstatus
        # hit_list = [line.split() for line in result if line]
        print(f'# Found {len(hit_list)} BLAST hits of the reference '
              f'pre-miRNA in the query')
        genes = pyfaidx.Fasta(query)
        # collect heuristic sequences
        # print('collecting sequences')
        # set border regions to include in analysis
        extraregion = 1000
        diff_dict = {}
        hit_at_start = False
        hit_at_end = False
        heuristic_fa = '{0}/candidate_region_{1}.out'.format(out, mirna_id)
        with open(heuristic_fa, 'w') as of:
            for hit in hit_list:
                chrom, start, end, strand, length = hit
                strand = strand.replace('plus', '+')
                strand = strand.replace('minus', '-')
                if strand == '-':
                    start, end = end, start
                header = ">{}|{}|{}|{}\n".format(chrom, start, end, strand)
                # print(header)
                start = int(start)
                end = int(end)
                if start > extraregion:
                    n_start = start - extraregion
                else:
                    # hit starts at the beginning of the chromosome
                    n_start = 0
                    hit_at_start = True
                n_end = end + extraregion
                if n_end > genes[chrom][-1].end:
                    # hit is at the end of the chromosome
                    n_end = genes[chrom][-1].end
                    hit_at_end = True
                if strand == '+':
                    sequence = genes[chrom][n_start:n_end].seq
                elif strand == '-':
                    sequence = genes[chrom][n_start:n_end].reverse.complement.seq
                header = ">{}|{}|{}|{}|{}|{}\n".format(chrom, start, end, strand, hit_at_start, hit_at_end)
                of.write(header)
                of.write(f'{sequence}\n')

        cm_results = []
        return_data = []
        cms_command = (
            'cmsearch --cpu {0} --noali -T {3} --incT {3} {1}/{2}.cm {4}'
            .format(cpu, models, mirna_id, cut_off, heuristic_fa)
        )
        cms_call = sp.Popen(
            cms_command, shell=True, stdout=sp.PIPE, stderr=sp.PIPE, encoding='utf8'
        )
        res, err = cms_call.communicate(mirna.pre)
        # delete temporary file
        if cleanup:
            os.remove(heuristic_fa)
        # read results
        if err:
            print(f'ERROR: {err}')
            sys.exit()
        for line in res.split('\n'):
            if line.startswith('  ('):
                if 'No hits detected that satisfy reporting thresholds' in line:
                    exitstatus = 'CMsearch found not hits'
                    return cm_results, exitstatus
                else:
                    data = line.strip().split()
                    cm_results.append(data)
        if not cm_results:
            exitstatus = 'CMsearch found no hits'
            return cm_results, exitstatus

        for data in cm_results:
            # print(data)
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
                print('Unexpected hit of CMsearch')
                continue
            if hit_at_start == 'True':
                print('# cmsearch hit for {} at the beginning of a chromosome in the query species'.format(mirna_id))
                hit_start = cm_start
                hit_end = cm_end
            elif hit_at_end == 'True':
                print('# cmsearch hit for {} at the end of a chromosome in the query species'.format(
                    mirna_id))
                hit_start = blast_start + (cm_start - extraregion) - 1
                hit_end = hit_start + (cm_end - cm_start) - 1
            else:
                hit_start = blast_start + (cm_start - extraregion) - 1
                hit_end = hit_start + (cm_end - extraregion) - 1
            # fill output variable
            return_data.append([blast_chrom, hit_start, hit_end, blast_strand, float(data[3])])
        if not return_data:
            exitstatus = 'No hits left after parsing of CMsearch results'
            return False, exitstatus

        cm_results, exitstatus = cmsearch_parser(return_data, cut_off, len_cut, mirna_id)

        if not cleanup:
            heur_results = '{0}/cmsearch_{1}.out'.format(out, mirna_id)
            with open(heur_results, 'w') as of:
                of.write('# CMsearch hits in query genome\n')
                for hit in cm_results:
                    out_d = [str(entry) for entry in cm_results[hit]]
                    of.write('\t'.join(out_d))
                    of.write('\n')

    # non heuristic mode (searches whole query genome with CM)
    else:
        cm_res_list = []
        cms_output = '{0}/cmsearch_{1}.out'.format(out, mirna_id)
        if not os.path.isfile(cms_output):
            # heuristic_cms = '{0}/cmsearch_heuristic_{1}.out'.format(out, mirna_id)
            print('# Running covariance model search for {}'.format(mirna_id))
            cms_command = (
                'cmsearch -T {5} --incT {5} --cpu {0} --noali '
                '-o {1} {2}/{3}.cm {4}'
                .format(cpu, cms_output, models, mirna_id, query, cut_off)
            )
            sp.run(cms_command, shell=True)

        else:
            print('# Found cm_search results at: {}. Using those'.format(cms_output))
        data = []
        with open(cms_output, 'r') as inf:
            for line in inf:
                if line.startswith('  ('):
                    if 'No hits detected that satisfy reporting thresholds' in line:
                        exitstatus = 'CMsearch found not hits'
                        return cm_results, exitstatus
                    data = line.strip().split()
                    h_chrom = data[5]
                    h_start = int(data[6])
                    h_end = int(data[7])
                    h_score = float(data[3])
                    h_strand = data[8]
                    cm_res_list.append([h_chrom, h_start, h_end, h_strand, h_score])
        if not data:
            exitstatus = 'cmSearch did not find anything'
            return cm_results, exitstatus

        cm_results, exitstatus = cmsearch_parser(cm_res_list, cut_off, len_cut, mirna_id)
    return cm_results, exitstatus


# cmsearch_parser: Parse the output of cmsearch while eliminating
#                  duplicates and filtering entries according to the
#                  defined cutoff.
def cmsearch_parser(cms, cmc, lc, mirid):
    """
    Parse the output of cmsearch while eliminating duplicates and filtering
    entries according to the defined cutoff.

    Parameters
    ----------
    cms     :   List with CMsearch hits [chrom, start, end, strand, score]
    cmc     :   Cutoff to decide which candidate hits should be included for the reverse BLAST search
    lc      :   Length cutoff
    mirid   :   miRNA ID

    Returns
    -------
    hits_dict : Dictionary containing hits of the cmsearch

    """
    # Output
    hits_dict = {}
    # Required for finding duplicates, stores hits per chromosome
    chromo_dict = {}
    cut_off = cmc

    for candidate_nr, hit in enumerate(cms, 1):
        h_chrom, h_start, h_end, h_strand, h_score = hit
        # blastparser expects the start to be the smaller number, will extract reverse complement if on - strand
        if h_start > h_end:
            start_tmp = h_start
            h_start = h_end
            h_end = start_tmp
        data = (
            '{0}_c{1}'.format(mirid, candidate_nr),
            h_chrom, h_start, h_end, h_strand, h_score
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
    exitstatus = 'Sucess'
    return hits_dict, exitstatus
