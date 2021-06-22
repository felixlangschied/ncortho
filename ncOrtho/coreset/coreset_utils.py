import sys
import re

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
                    and line.split('\t')[2] == 'gene'
                    and 'gene_biotype=protein_coding' in line.split()[-1]
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
