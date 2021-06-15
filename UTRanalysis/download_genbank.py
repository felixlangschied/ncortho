from Bio import Entrez
import xml.etree.ElementTree as ET
import wget
from ftplib import FTP
import os
import json


# define I/O
phylprofile_results = '/home/felixl/project/ncOrtho/analyses/diff_sets/hsa_chromonly_heuristic_dustfiltered/rodent_missing/fDOG/ingos_run.phyloprofile'
# phylprofile_results = r'ingos_run.phyloprofile'
outputdir = '/share/project/felixl/ncOrtho/data/UTR/genbank_files'
# outputdir = '.'
# cats_to_extract = ['Primates', 'Rodentia']
cats_to_extract = ['Primates', 'Rodentia', 'Artiodactyla']

# get scientific name from taxid
def name_from_taxid(taxid, tax_units, email):
    return_names = {}
    Entrez.email = email
    handle = Entrez.efetch(db="taxonomy", id=taxid, retmode='xml')
    raw_data = handle.read().decode()
    root = ET.fromstring(raw_data)
    handle.close()
    for c, taxon in enumerate(root):
        for child in taxon:
            if child.tag == 'TaxId':
                this_taxid = child.text
            if child.tag == 'ScientificName':
                name = child.text
                name = name.replace(' ', '_')
            if child.tag == 'Lineage':
                isIngroup = []
                for group in tax_units:
                    if group in child.text:
                        isIngroup.append(True)
                    else:
                        isIngroup.append(False)
                if any(isIngroup):
                    return_names[this_taxid] = name
    return return_names

# download rna.genbank files from Refseq FTP
def extract_GCF(tax_name, output):
    ftp = FTP('ftp.ncbi.nlm.nih.gov')
    ftp.login()
    try:
        ftp.cwd("/genomes/refseq/vertebrate_mammalian/{}/latest_assembly_versions/".format(tax_name))
        assemblies = ftp.nlst()
        assembly = assemblies[0]

        download_url = (
            'ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/vertebrate_mammalian/{0}/latest_assembly_versions/{1}/{1}_rna.gbff.gz'
                .format(tax_name, assembly)
        )
    except:
        print('\n# {} not found on RefSeq FTP'.format(tax_name))
        new_name = '_'.join(tax_name.split('_')[0:2])
        print('# Trying with {}'.format(new_name))
        try:
            ftp.cwd("/genomes/refseq/vertebrate_mammalian/{}/latest_assembly_versions/".format(new_name))
            assemblies = ftp.nlst()
            assembly = assemblies[0]
            download_url = (
                'ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/vertebrate_mammalian/{0}/latest_assembly_versions/{1}/{1}_rna.gbff.gz'
                    .format(new_name, assembly)
            )
        except:
            print('# Still could not find {} on RefSeq FTP. Skipping..'.format(new_name))
            return None

    os.chdir(output)
    if not os.path.isfile(f'{output}/{assembly}_rna.gbff') and not os.path.isfile(f'{output}/{assembly}_rna.gbff.gz'):
        print('\n# Starting download of rna.genbank file for {}'.format(tax_name))
        wget.download(download_url)
    else:
        print('\n# rna.genbank file already existant. Skipping Download..')
    return assembly

# main
# prepare phyloprofile output
taxids = set()
with open(phylprofile_results, 'r') as fh:
    # skip header
    next(fh)
    for line in fh:
        data = line.strip().split()
        taxid = data[1].replace('ncbi', '')
        taxids.add(taxid)
taxids = ', '.join(list(taxids))

# parse Taxonomy ids
tax_names = name_from_taxid(taxids, cats_to_extract, 'langschied@bio.uni-frankfurt.de')
# with open(f'{outputdir}/taxname_dict.json', 'w') as of:
#     json.dump(tax_names, of)

# start download
with open(f'{outputdir}/taxid_name_assembly.tsv', 'w') as of:
    for taxid in tax_names:
        name = tax_names[taxid]
        assembly = extract_GCF(name, outputdir)
        of.write('{}\t{}\t{}\n'.format(taxid, name, assembly))



