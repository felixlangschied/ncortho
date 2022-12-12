import sys
import re
import yaml


def vprint(s, verbose):
    if verbose:
        print(s, flush=True)

# Parse an Ensembl GTF file to store the coordinates for each protein-coding gene in a
# dictionary
def gtf_parser(gtf):
    species_name = gtf.split('/')[-1].split('.')[0]
    chr_dict = {}
    chromo = ''

    # with open(inpath) as infile, open(outpath, 'wb') as outfile:
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
                    raise ValueError('Unknown identifier type "{}"'.format(id_type))

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


def table_parser(path):
    chromo = ''
    chr_dict = {}
    with open(path) as fh:
        for line in fh:
            if line.startswith('#'):
                continue
            contig, start, end, strand, geneid = line.strip().split('\t')
            if contig != chromo:
                i = 1
                chromo = contig
            chr_dict[geneid] = (contig, i)
            try:
                chr_dict[contig][i] = (geneid, start, end, strand)
            except KeyError:
                chr_dict[contig] = {i: (geneid, start, end, strand)}
            i += 1
    return chr_dict


def parse_annotation(path, idtyp):
    ft = path.split('.')[-1]
    if ft == 'gtf':
        rd = gtf_parser(path)
    elif ft in ['gff3', 'gff']:
        rd = gff_parser(path, idtyp)
    elif ft in ['tsv', 'txt']:
        rd = table_parser(path)
    else:
        raise ValueError(f'File type "{idtyp}" not valid as reference annotation. Expecting .gff3, .gff or .gtf')
    # print('Done')
    return rd


###############################################################################

def parse_yaml(path):
    p_dict = {}  # {<name>: {'orthologs': <path>, 'genome': <path>, 'annotation': <path>}}
    paths = []
    refdict = {}
    with open(path, 'r') as param_handle:
        params = yaml.load_all(param_handle, Loader=yaml.FullLoader)
        for entry in params:
            in_type = entry.pop('type')
            if in_type == 'reference':
                refdict['annotation'] = entry['annotation']
                refdict['genome'] = entry['genome']
            elif in_type == 'core':
                species = entry.pop('name')
                p_dict[species] = entry
            paths.extend([entry['annotation'], entry['genome']])
        return p_dict, refdict, paths


def read_pairwise_orthologs(pathdict):
    od = {}
    for taxon, pd in pathdict.items():
        taxd = {}
        orthofile = pd['orthologs']
        with open(orthofile) as omafile:
            for line in omafile:
                if line.startswith('#'):
                    continue
                ref, core = line.strip().split()[:2]
                if not ref in taxd:
                    taxd[ref] = []
                taxd[ref].append(core)
        od[taxon] = taxd
    return od


def read_mirnas(path):
    with open(path) as mirfile:
        mirnas = [
            line.split() for line in mirfile.readlines()
            if not line.startswith('#')
        ]
    # print('Done')
    return mirnas





def mirna_position(mirlist):
    mirid, rawchrom, start, end, strand = mirlist[:5]
    if 'chr' in rawchrom:
        chromo = rawchrom.split('chr')[1]
    else:
        chromo = rawchrom
    return mirid, chromo, int(start), int(end), strand





