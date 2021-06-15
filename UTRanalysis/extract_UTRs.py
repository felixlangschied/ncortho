from Bio import SeqIO
import json
import os


# tax_path = '/share/project/felixl/ncOrtho/data/UTR/taxid_name_assembly.tsv'
tax_path = r'C:\Users\felix\project\UTRanalysis\taxid_name_assembly.tsv'
genbank_files_dir = '/share/project/felixl/ncOrtho/data/UTR/genbank_files'
out_dir = '/share/project/felixl/ncOrtho/data/UTR/utr_jsons'


with open(tax_path, 'r') as fh:
    for line in fh:
        taxid, name, assembly = line.strip().split('\t')
        outpath = '{}/{}_UTRs.json'.format(out_dir, assembly)
        if assembly == 'None':
            continue
        elif os.path.isfile(outpath):
            continue
        genbank_loc = '{}/{}_rna.gbff'.format(genbank_files_dir, assembly)

        out_dict = {}
        with open(genbank_loc) as handle:
            for record in SeqIO.parse(handle, "genbank"):
                nuc_seq = record.seq
                for feature in record.features:
                    if feature.type =='CDS':
                        prot_id = feature.qualifiers['protein_id'][0]
                        loc = feature.location
                        prot_dict = {}
                        prot_dict['five_UTR'] = str(nuc_seq[:loc.start] )
                        prot_dict['three_UTR'] = str(nuc_seq[loc.end:])
                        out_dict[prot_id] = prot_dict
        print('# UTR extraction for {} done'.format(name))


        if not os.path.isdir(out_dir):
            os.mkdir(out_dir)
        with open(outpath, 'w') as of:
            json.dump(out_dict, of)