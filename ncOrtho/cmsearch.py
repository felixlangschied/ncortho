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

import glob
import pyfaidx
import subprocess as sp
import sys
import os


def cmsearcher(mirna, cm_cutoff, cpu, msl, models, query, blastdb, out, cleanup, heuristic):
    """
    Parameters
    ----------
    mirna       :   Mirna object
    cm_cutoff   :   Minimum bit score of the CMsearch hit relative to the
                    bit score that results from searching in the reference genome with the same model
    cpu         :   Number of cores to use
    msl             Minimum length of CMsearch hit relative to the reference pre-miRNA length
    models      :   Path to the directory that contains the CMs
    query       :   Path to the query genome
    qblast      :   Optional Path to a BLASTdb of the query genome
    out         :   Path to the output directory
    cleanup     :   Delete intermediate files True/False
    heuristic   :   Should heuristic mode be used and if yes, set parameters (True/False, evalue, minlength)

    Returns
    -------
    cm_results : Dictionary containing hits of the cmsearch

    """
    mirna_id = mirna.name

    # blastdb = query.replace('.fa', '')
    # Calculate the bit score cutoff.
    cut_off = mirna.bit * cm_cutoff
    # Calculate the length cutoff.
    len_cut = len(mirna.pre) * msl
    blast_len_cut = len(mirna.pre) * heuristic[2]
    # if no model exists for mirna return False
    if not os.path.isfile('{}/{}.cm'.format(models, mirna_id)):
        print('# No model found for {}. Skipping..'.format(mirna_id))
        return False
    # Perform covariance model search.
    # Report and inclusion thresholds set according to cutoff.

    if heuristic[0]:
        # extract candidate regions
        print('# Identifying candidate regions for cmsearch heuristic')
        tmp_out = '{0}/tmp_blast_{1}.out'.format(out, mirna_id)
        with open(tmp_out, 'w') as inf:
            inf.write(">{}\n{}\n".format(mirna_id, mirna.pre))
        blast_command = (
            'blastn -evalue {3} -task blastn -db {0} -query {1} '
            '-num_threads {2} -outfmt "6 sseqid sstart send sstrand length"'.format(blastdb, tmp_out, cpu, heuristic[1])
        )
        # print(blast_command)
        res = sp.run(blast_command, shell=True, capture_output=True)
        os.remove(tmp_out)
        if res.stdout:
            print('Blast step finished')
        else:
            cm_results = False
            return cm_results
        result = res.stdout.decode('utf-8').split('\n')
        hit_list = [line.split() for line in result if line and float(line.split()[-1]) >= blast_len_cut]
        # hit_list = [line.split() for line in result if line]
        print(f'# Found {len(hit_list)} BLAST hits of the reference '
              f'pre-miRNA in the query')
        genes = pyfaidx.Fasta(query)
        # collect heuristic sequences
        print('collecting sequences')
        # set border regions to include in analysis
        extraregion = 1000
        diff_dict = {}
        hit_at_start = False
        hit_at_end = False
        heuristic_fa = '{0}/heuristic_{1}.out'.format(out, mirna_id)
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
                    hit_at_start, hit_at_end = data[0].split('|')[4:6]
                    blast_strand = data[0].split('|')[3]
                    cm_start, cm_end = map(int, data[7:9])
                    # cm_start, cm_end = map(int, data[5:7])
                    cm_strand = data[9]
                    if cm_strand == '-':
                        # reverse hits should not be possible since blasthits were extracted
                        # as reverse complement if they were on the minus strand
                        continue
                    # if not cm_strand == blast_strand:
                    #     print('Rejecting cmsearch hit. Not on same strand as BLAST hit')
                    #     break
                    # print(type(hit_at_start), hit_at_end)
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
                    # hit_start = blast_end - cm_end
                    # hit_end = blast_end - cm_start
                    data[0] = blast_chrom
                    data[7] = hit_start
                    data[8] = hit_end
                    data[9] = blast_strand
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
def cmsearch_parser(cms, cmc, lc, mirid):
    """
    Parse the output of cmsearch while eliminating duplicates and filtering
    entries according to the defined cutoff.

    Parameters
    ----------
    cms     :   Path to cmsearch output
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
                h_chrom = hit[0]
                h_start = hit[7]
                h_end = hit[8]
                h_strand = hit[9]
                h_score = hit[14]
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
    return hits_dict
