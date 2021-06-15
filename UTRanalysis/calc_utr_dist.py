import json
import subprocess as sp
from Bio import pairwise2
from Bio.Phylo.TreeConstruction import DistanceCalculator
from Bio import AlignIO
import pandas as pd
import numpy as np
import os

# I/O
# ref_path = '/share/project/felixl/ncOrtho/data/UTR/utr_jsons/GCF_000001405.39_GRCh38.p13_UTRs.json'
json_path = '/share/project/felixl/ncOrtho/data/UTR/utr_jsons'
phylo_output = '/home/felixl/project/ncOrtho/analyses/diff_sets/hsa_chromonly_heuristic_dustfiltered/rodent_missing/fDOG/ingos_run.phyloprofile'
info_path = '/share/project/felixl/ncOrtho/data/UTR/taxid_name_assembly.tsv'
outdir = '/share/project/felixl/ncOrtho/data/UTR/output'

# set this to True if you want to recalculate all alignments
force_aln = False

# protein_ids = ['O75896']
ref_taxa = '9606'
# query_taxa = 10090
# fastme dna models. Alternatively, put 'biopython' as method to use the DistanceCalculator
fastme_models = ['p-distance', 'RY symmetric', 'RY', 'JC69', 'K2P', 'F81', 'F84', 'TN93', 'LogDet']
# method = 'F84'
method = 'biopython'

def utr_extract(utr_dict_path, ref_protid, taxid, ensembl_2_refseq, UTR='3'):
    with open(utr_dict_path, 'r') as fh:
        utr_dict = json.load(fh)
    try:
        refseq_id = ensembl_2_refseq[str(taxid)][ref_protid]
    except KeyError:
        print('# fDOG found no ortholog of {} in species with taxid {}'.format(ref_protid, taxid))
        return None
    if UTR == '3':
        try:
            utr_seq = utr_dict[refseq_id]['three_UTR']
            return utr_seq
        except KeyError:
            print('# Did not find UTR for protein {} in species with taxid {}'.format(refseq_id, taxid))
            return None
    elif UTR == '5':
        try:
            utr_seq = utr_dict[refseq_id]['five_UTR']
            return utr_seq
        except KeyError:
            print('# Did not find UTR for protein {} in species with taxid {}'.format(refseq_id, taxid))
            return None


def parse_phyloprofile(phylo_path):
    # dictionary used for mapping of ids between refseq and ensembl
    map_dict = {}
    # dictionary to map taxid to taxname in the format of: {'10036': 'MESAU@10036@1'}
    # taxid_2_taxname = {}
    # collect all protein_ids
    prot_ids = set()
    with open(phylo_path, 'r') as infile:
        # skip header
        next(infile)
        for line in infile:
            ensembl_id, taxa, refseq_id, isRefOrth = line.strip().split()[2].split('|')
            # skip non reference orthologs
            # TODO: calculate UTR dist for all homologs
            if not isRefOrth == '1':
                continue
            prot_ids.add(ensembl_id)
            taxid = taxa.split('@')[1]
            tax_dict = {ensembl_id: refseq_id}
            if taxid not in map_dict:
                map_dict[taxid] = tax_dict
            else:
                map_dict[taxid].update(tax_dict)
    return map_dict, list(prot_ids)


def calc_distmat(dna_mod, fasta_path, faln=False):
    # align UTRs using MUSCLE
    print('# Starting alignemnt')
    aln_out = fasta_path.replace('.fa', '.phys')
    aln_cmd = (
        f'muscle -phyi -in {fasta_path} -out {aln_out}'
    )
    if not os.path.isfile(aln_out) or faln:
        sp.run(aln_cmd, shell=True)
    else:
        print('# Found alignment at {}. Using this one instead..'.format(aln_out))
    matrix_out = aln_out.replace('.phys', '.matrix')
    # calculate distance matrix
    if dna_mod in fastme_models:
        # run fastME
        # change -d for dna model used. K2P for kimura. default: F84
        print('# Calculating distance of UTRs using the DNA model: {}'.format(dna_mod))
        fastme_cmd = f'fastme -c -i {aln_out} -O {matrix_out} -d {dna_mod}'
        sp.run(fastme_cmd, shell=True)
        return matrix_out
    elif dna_mod == 'biopython':
        aln = AlignIO.read(open(f'{aln_out}'), 'phylip')
        calculator = DistanceCalculator('identity')
        dm = calculator.get_distance(aln)
        return dm

# TODO: cannot handle NA cases atm
def parse_matrix(mat_path, ref_tax):
    tax_list = []
    prot_id = mat_path.split(os.sep)[-1].split('.')[0]
    with open(mat_path, 'r') as fh:
        mat_len = int(next(fh))
        data_ar = np.empty([mat_len, mat_len])
        c = 0
        for line in fh:
                data = line.strip().split()
                if data:
                    tax_list.append(data[0])
                    data_ar[c, :] = data[1:]
                    c += 1
    df = pd.DataFrame(data_ar, columns=tax_list, index=tax_list)
    ref_series = df[f'ncbi{ref_tax}'].rename(prot_id)
    return ref_series


# create output directory
if not os.path.isdir(outdir):
    os.mkdir(outdir)

# read phyloprofile output
ensembl_2_refseq, reference_protids = parse_phyloprofile(phylo_output)
print('# Finisthed parsing the PhyloProfile output')


# TODO: remove this break
# protein_ids = reference_protids[:2]
# reference_protids = ['P52926']
# exit()

# loop start
single_pylo_outs = []
taxids_sampled = set()
for prot_id in reference_protids:
    print('### Starting analysis of {} ###'.format(prot_id))
    # create output dir per protein:
    prot_out = '{}/{}'.format(outdir, prot_id)
    if not os.path.isdir(prot_out):
        os.mkdir(prot_out)
    # create multifasta of the UTRs of one protein
    fasta_out = '{}/{}.fa'.format(prot_out, prot_id)
    with open(fasta_out, 'w') as of, open(info_path, 'r') as fh:
        for line in fh:
            taxid, name, assembly = line.strip().split()
            path_2_json = '{}/{}_UTRs.json'.format(json_path, assembly)
            if assembly == 'None':
                continue
            elif not os.path.isfile(path_2_json):
                print('# No UTR json found for {}. Skipping..'.format(name))
                continue
            else:
                taxids_sampled.add(taxid)
                utr = utr_extract(path_2_json, prot_id, taxid, ensembl_2_refseq)
                if utr and utr is not None:
                    of.write('>ncbi{}\n'.format(taxid))
                    of.write(utr)
                    of.write('\n')
    print('# Crated multifasta file for {}'.format(prot_id))

    # calculate distance matrix
    dm_out = calc_distmat(method, fasta_out, force_aln)
    print('# Finished distance calculation for {}\n'.format(prot_id))
    # print('# Writing UTR PhyloProfile file')
    if method == 'biopython':
        pp_out = '{}/utr_{}.phyloprofile'.format(prot_out, prot_id)
        single_pylo_outs.append(pp_out)
        with open(phylo_output, 'r') as infile, open(pp_out, 'w') as of:
            # skip header
            header = next(infile).split()
            # write header
            header[-2] = '3_UTR'
            header[-1] = 'TBD'
            header = '\t'.join(header)
            of.write(header)
            of.write('\n')
            for line in infile:
                ref_prot, taxid, info, val1, val2 = line.strip().split()
                if info.split('|')[-1] != '1':
                    continue
                elif taxid.replace('ncbi', '') in taxids_sampled:
                    if ref_prot == prot_id:
                        try:
                            utr_score = 1 - dm_out[f'ncbi{ref_taxa}', taxid]
                            outstring = '{}\t{}\t{}\t{}\t0\n'.format(ref_prot, taxid, info, utr_score)
                            of.write(outstring)
                        except ValueError:
                            outstring = '{}\t{}\t{}\t0\t0\n'.format(ref_prot, taxid, info)
                            of.write(outstring)
    # extract distances to the reference species UTRs and concat them to dataframe
    # elif method in fastme_models:
    #     df = []
    #     for prot_id in reference_protids:
    #         prot_mat_path = '{0}{2}{1}{2}{1}.matrix'.format(outdir, prot_id, os.sep)
    #         series = parse_matrix(prot_mat_path, taxid)
    #         if not df:
    #             df = pd.DataFrame(series)
    #         else:
    #             pd.concat([df, series])
        # df.to_csv(f'{outdir}/dist_df.csv')

# write total output
total_out = '{}/utr_dist.phyloprofile'.format(outdir)
first = True
with open(total_out, 'w') as of:
    for pp_f in single_pylo_outs:
        with open(pp_f, 'r') as fh:
            header = next(fh)
            if first:
                of.write(header)
            else:
                of.write(fh.read())
            first = False
print('# Analysis finished. Results are located at:')
print(total_out)










