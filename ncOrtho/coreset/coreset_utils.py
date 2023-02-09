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
            if line.startswith('#') or not 'gene_biotype=protein_coding' in line:
                continue
            linedata = line.strip().split('\t')
            if linedata[2] != 'gene':
                continue

            if id_type == 'GeneID':
                geneid = re.split('[;,]', linedata[-1].split(f'{id_type}:')[1])[0]
            elif id_type in ['ID', 'Name', 'gene_id', 'CDS']:
                geneid = linedata[-1].split(f'{id_type}=')[1].split(';')[0]
                if id_type in ['ID', 'CDS']:
                    geneid = '-'.join(geneid.split('-')[1:])
            else:
                raise ValueError('Unknown identifier type "{}"'.format(id_type))

            contig = linedata[0]
            start = int(linedata[3])
            end = int(linedata[4])
            strand = linedata[6]
            if contig != chromo:  # reset counter
                i = 1
                chromo = contig
            chr_dict[geneid] = (contig, i)
            if contig not in chr_dict:
                chr_dict[contig] = {}
            chr_dict[contig][i] = (geneid, start, end, strand)

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
                chr_dict[contig][i] = (geneid, int(start), int(end), strand)
            except KeyError:
                chr_dict[contig] = {i: (geneid, int(start), int(end), strand)}
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


def gene_from_cds(gff):
    """
    returns:    d = {'WP_23898234': 'CDK'}
    """
    # cds2parent = {}
    cds_mrna_map = {}

    cds2gene = {}
    orphan = []
    with open(gff) as fh:
        for line in fh:
            if line.startswith('#'):
                continue
            data = line.strip().split('\t')
            if data[2] != 'CDS':
                continue
            cdsid = data[-1].split(';')[0].split('.')[0].replace('ID=cds-', '').replace('ID=id-', '')

            parent = data[-1].split('Parent=')[1].split(';')[0].split('.')[0]
            if parent.startswith('gene-') or parent.startswith('id-'):
                cds2gene[cdsid] = parent.split(';')[0].replace('gene-', '').replace('id-', '')
            else:
                cds_mrna_map[cdsid] = parent
                cds_mrna_map[parent] = cdsid

    with open(gff) as fh:
        for line in fh:
            if line.startswith('#'):
                continue
            data = line.strip().split('\t')
            if data[2] != 'mRNA':
                continue
            mrnaid = data[-1].split('ID=')[1].split('.')[0]
            parent = data[-1].split('Parent=')[1].split('.')[0]
            if mrnaid not in cds_mrna_map:
                continue
            cdsid = cds_mrna_map[mrnaid]
            if parent.startswith('gene-'):
                cds2gene[cdsid] = parent.split(';')[0].replace('gene-', '')
            else:
                orphan.append(parent)
    return cds2gene


def pairwise_orthologs_from_cds(cd, refanno):
    """
    returns    d = {corespec: {refprot: [ortho1, ortho2, ...], ...}, ...}
    """
    ref_cds2gene = gene_from_cds(refanno)
    od = {}
    for cspec in cd:
        cds2gene = gene_from_cds(cd[cspec]['annotation'])

        od[cspec] = {}
        orthopath = cd[cspec]['orthologs']
        with open(orthopath) as fh:
            for line in fh:
                cds_refprot, cds_orthoprot = line.strip().split()
                try:
                    refprot = ref_cds2gene[cds_refprot]
                except KeyError as e:
                    # print('reference', e)
                    # IDs of "longest proteins" that were removed in later versions of the RefSeq annotation.
                    # Discrepency possibly due to versioning problems
                    continue
                try:
                    orthoprot = cds2gene[cds_orthoprot]
                except KeyError as e:
                    print(cspec, e)
                    continue

                if not refprot in od[cspec]:
                    od[cspec][refprot] = []
                od[cspec][refprot].append(orthoprot)
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





