"""
Created on Fri Nov 20 14:24:19 2020
@author: felixl

Download a set of genomes from EnsemblDB or NCBI
based on NCBI accession number or Ensembl genus/species name
e.g a file containing a list of accession numbers like:
GCA_000001635.4
GCF_000001405.39

or a file containing Ensembl genus/species names like:
Mus_musculus
Homo_sapiens

USAGE:
    flavorSavor -i /path/to/IDs2download -o /path/to/output -f flavor

use genome_downloader.py -h for a list of all available flavors
"""

from ftplib import FTP
import re
import wget
import os
import subprocess as sp
import sys
import argparse
import textwrap


# download data from ensembl FTP based on ensembl ids
def extract_ensembl(ids, output, raw_flavor, naming, dir_mode):
    ids = [eid.lower() for eid in ids]
    ftp = FTP('ftp.ensembl.org')
    ftp.login()
    ftp.cwd("/pub/current_fasta")
    root_dirs = ftp.nlst()
    flavor = raw_flavor.replace('ensembl_', '')
    count = -1

    for id in ids:
        count += 1
        if id in root_dirs:
            if flavor in ('dna', 'pep', 'cds', 'dna_sm'):
                direct = 'current_fasta'
                if flavor == 'dna':
                    ftp.cwd('/pub/{}/{}/{}'.format(direct, id, flavor))
                    files = ftp.nlst()
                    r = re.compile(".*dna.toplevel.fa.gz")
                    target = list(filter(r.match, files))[0]
                elif flavor == 'dna_sm':
                    flavor = 'dna'
                    ftp.cwd('/pub/{}/{}/{}'.format(direct, id, flavor))
                    files = ftp.nlst()
                    r = re.compile(".*dna_sm.toplevel.fa.gz")
                    target = list(filter(r.match, files))[0]
                elif flavor == 'pep':
                    ftp.cwd('/pub/{}/{}/{}'.format(direct, id, flavor))
                    files = ftp.nlst()
                    r = re.compile('.*pep.all.fa.gz')
                    target = list(filter(r.match, files))[0]
                elif flavor == 'cds':
                    ftp.cwd('/pub/{}/{}/{}'.format(direct, id, flavor))
                    files = ftp.nlst()
                    r = re.compile('.*cds.all.fa.gz')
                    target = list(filter(r.match, files))[0]
                try:
                    t_location = (
                        'ftp://ftp.ensembl.org/pub/{0}/{1}/{2}/{3}'
                        .format(direct, id, flavor, target)
                    )
                    os.chdir(output)
                    print('\n# Starting download for {}'.format(id))
                    wget.download(t_location)
                    # rename files if option selected
                    if naming:
                        name = naming[count]
                        if dir_mode:
                            if not os.path.isdir('{}/{}'.format(output, name)):
                                cmd = 'mkdir {}'.format(name)
                                sp.run(cmd, shell=True)
                            file_name = '{}.fa.gz'.format(name)
                            cmd = 'mv {} {}/{}'.format(target, name, file_name)
                            sp.run(cmd, shell=True)
                        else:
                            file_name = '{}.fa.gz'.format(name)
                            cmd = 'mv {} {}'.format(target, file_name)
                            sp.run(cmd, shell=True)
                    elif not naming and dir_mode:
                        print('\n # No column with names detected in the input file. Unable to run in dir_mode')
                except:
                    print('\n# {} not found'.format(target))
            elif flavor in ('gtf', 'gff3'):
                direct = 'current_{}'.format(flavor)
                #print('/pub/{}/{}'.format(direct, id))
                ftp.cwd('/pub/{}/{}'.format(direct, id))
                files = ftp.nlst()
                for file in files:
                    if (
                            '.gz' in file
                            and not ('abinitio.gtf.gz' and 'chr.gtf.gz') in file
                    ):
                        target = file
                try:
                    t_location = (
                        'ftp://ftp.ensembl.org/pub/{0}/{1}/{2}'
                            .format(direct, id, target)
                    )
                    os.chdir(output)
                    print('\n# Starting download for {}'.format(id))
                    wget.download(t_location)
                    # rename files if option selected
                    if naming:
                        name = naming[count]
                        if dir_mode:
                            if not os.path.isdir('{}/{}'.format(output, name)):
                                cmd = 'mkdir {}'.format(name)
                                sp.run(cmd, shell=True)
                            file_name = '{}.gtf.gz'.format(name)
                            cmd = 'mv {} {}/{}'.format(target, name, file_name)
                            sp.run(cmd, shell=True)
                        else:
                            file_name = '{}.gtf.gz'.format(name)
                            cmd = 'mv {} {}'.format(target, file_name)
                            sp.run(cmd, shell=True)
                    elif not naming and dir_mode:
                        print('\n# No column with names detected in the input file. Unable to run in dir_mode')
                except:
                    print('\n# {} not found'.format(target))
        else:
            print('\n# {} not found on the Ensembl FTP-server'.format(id))


# download genomes from NCBI based on GCF/GCA id
def extract_GCF(ids, output, flavor, naming, dir_mode):
    ftp = FTP('ftp.ncbi.nlm.nih.gov')
    ftp.login()
    ftp.cwd("/genomes/all/")
    count = -1

    for id in ids:
        count += 1
        ftp.cwd('/genomes/all/')
        id = id.split('.')[0]
        typ = id.split('_')[0]
        id = id.split('_')[1]
        part1 = id[0:3]
        part2 = id[3:6]
        part3 = id[6:9]

        url = '/'.join([typ, part1, part2, part3])
        ftp.cwd(url)
        assemblies = ftp.nlst()
        # will download the top file
        target = assemblies[0]
        if len(assemblies) > 1:
            print('# Found multiple assemblies for {}.\n# Downloading most recent: {}'
                  .format(id, assemblies[0])
                  )
        ftp.cwd(target)
        files = ftp.nlst()

        # find file types that can be downloaded (flavors)
        # ftypes = [file.split('.')[-3:] for file in files if '.gz' in file]
        # ftypes = ['.'.join(parts) for parts in ftypes]
        # ftypes = [ftype.split('_')[2:] for ftype in ftypes]
        # ftypes = ['_'.join(parts) for parts in ftypes]
        # print(ftypes)
        # print('\n'.join(ftypes))

        r = re.compile('.*' + flavor)
        target_file = list(filter(r.match, files))[0]
        download_url = (
            'ftp://ftp.ncbi.nlm.nih.gov/genomes/all/{}/{}/{}'
                .format(url, target, target_file)
        )
        os.chdir(output)
        try:
            print('\n# Starting download of {} for: {}'.format(flavor, target))
            wget.download(download_url)
            if naming:
                name = naming[count]
                print(name)
                if dir_mode:
                    if not os.path.isdir('{}/{}'.format(output, name)):
                        cmd = 'mkdir {}'.format(name)
                        sp.run(cmd, shell=True)
                    file_name = '{}_{}'.format(name, flavor)
                    cmd = 'mv {}_{} {}/{}'.format(target, flavor, name, file_name)
                    sp.run(cmd, shell=True)
                else:
                    file_name = '{}_{}'.format(name, flavor)
                    cmd = 'mv {}_{} {}'.format(target, flavor, file_name)
                    sp.run(cmd, shell=True)
            elif not naming and dir_mode:
                print('\n# No column with names detected in the input file. Unable to run in dir_mode')

        except:
            print('\n# {} not found'.format(target_file))


# Allow boolean argument parsing
def str2bool(v):
    if isinstance(v, bool):
        return v
    if v.lower() in ('yes', 'true', 't', 'y', '1'):
        return True
    elif v.lower() in ('no', 'false', 'f', 'n', '0'):
        return False
    else:
        raise argparse.ArgumentTypeError('Boolean value expected.')


# RUN PROGRAM
def main():
    # noinspection PyTypeChecker
    parser = argparse.ArgumentParser(
        prog='python flavorSavor.py',
        description='Download a set of genomes from EnsemblDB or NCBI '
                    'based on NCBI accession number or Ensembl genus/species name ',
        formatter_class=argparse.RawTextHelpFormatter
    )
    parser.add_argument(
        '-i', '--input', metavar='<path>', type=str,
        help='Path to list of NCBI-Accessions (e.g. GCA_000001635.4) or ensembl genus/species names (e.g. Mus_musculus). '
             'Can have an optional tab delimited column that speciefies a name that the downloaded file should receive ',
    )
    parser.add_argument(
        '-o', '--output', metavar='<path>', type=str,
        help='Output path'
    )
    parser.add_argument(
        '-f', '--flavor', metavar='<str>', type=str,
        help=textwrap.dedent('''\
        Type of data you want to download.

        # List of available flavors when downloading from NCBI:
        genomic.fna.gz
        genomic.gbff.gz
        genomic.gff.gz
        protein.faa.gz
        protein.gpff.gz
        feature_table.txt.gz
        wgsmaster.gbff.gz
        cds_from_genomic.fna.gz
        rna_from_genomic.fna.gz
        feature_count.txt.gz
        translated_cds.faa.gz
        genomic.gtf.gz
        genomic_gaps.txt.gz

        # If you want to download from Ensembl use a flavor from this list:
        ensembl_dna
        ensembl_dna_sm
        ensembl_pep
        ensembl_gtf
        ensembl_gff3
        ensembl_cds

        # Currently not supported ensembl formats:
        ensembl_cdna
        ensembl_cds
        ensembl_dna_index
        ensembl_ncrna
        ...

        ''')
    )
    parser.add_argument(
        '-u', '--unpack', type=str2bool, metavar='True/False', nargs='?', const=True, default=True,
        help="Set to False if you don't want to unpack the downloaded files straight away (Default=True)"
    )
    parser.add_argument(
        '-d', '--dir_mode', type=str2bool, metavar='True/False', nargs='?', const=False, default=False,
        help="Downloads each file into a new directory. "
             "Only applicable with name column in the input file (Default=False)"
    )
    # Show help when no arguments are added.
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)
    else:
        args = parser.parse_args()
    # parse input
    output = os.path.abspath(args.output)
    input = args.input
    flavor = args.flavor
    unpack = args.unpack
    dir_mode = args.dir_mode

    # MAIN BODY
    # Check if output folder exists or create it otherwise
    if not os.path.isdir(output):
        print('# Creating output folder')
        cmd = 'mkdir {}'.format(output)
        try:
            sp.run(cmd, shell=True, check=True, stderr=sp.PIPE)
        except sp.CalledProcessError:
            print('# Could not create output folder at:\n'
                  '{}\n'
                  'Exiting..'.format(output))

    # Verify flavor
    ncbi_flavors = ['genomic.fna.gz', 'genomic.gbff.gz', 'genomic.gff.gz',
                    'protein.faa.gz', 'protein.gpff.gz',
                    'feature_table.txt.gz', 'wgsmaster.gbff.gz',
                    'cds_from_genomic.fna.gz', 'rna_from_genomic.fna.gz',
                    'feature_count.txt.gz', 'translated_cds.faa.gz',
                    'genomic.gtf.gz', 'genomic_gaps.txt.gz']

    ensembl_supported = ['ensembl_dna', 'ensembl_pep', 'ensembl_gtf', 'ensembl_gff3', 'ensembl_cds', 'ensembl_dna_sm']
    ensembl_unsupported = ['ensembl_cdna', 'ensembl_cds',
                           'ensembl_dna_index', 'ensembl_ncrna']

    if not flavor in ncbi_flavors + ensembl_supported:
        print('# Unknown flavor: {}'.format(flavor))
        sys.exit()
    elif flavor in ensembl_unsupported:
        print('# Flavor not yet supported: {}'.format(flavor))
        sys.exit()

    # Read input
    with open(input, 'r') as file:
        id_list = []
        name_list = []
        naming = []
        data = file.read().strip().split('\n')
        for line in data:
            id_list.append(line.split('\t')[0])
            if len(line.split('\t')) == 2:
                name_list.append(line.split('\t')[1])
        if name_list and len(name_list) == len(id_list):
            naming = name_list
        # convert spaces to underscores
        id_list = [id.replace(' ', '_') for id in id_list]
    # Download
    if 'ensembl_' in flavor:
        extract_ensembl(id_list, output, flavor, naming, dir_mode)
    else:
        extract_GCF(id_list, output, flavor, naming, dir_mode)
    # Unpack the downloaded files
    if unpack and not dir_mode:
        print('\n# Trying to unpack downloaded files..')
        os.chdir(output)
        cmd = 'gunzip *.gz'
        sp.run(cmd, shell=True)
    elif unpack and dir_mode:
        print('\n# Trying to unpack downloaded files..')
        if not naming:
            print('# Unable to run in dir_mode. Trying to unpack at the output directory ')
            os.chdir(output)
            cmd = 'gunzip *.gz'
            sp.run(cmd, shell=True, check=True, stderr=sp.PIPE)
        else:
            for name in naming:
                path = '{}/{}'.format(output, name)
                if os.path.isdir(path):
                    os.chdir(path)
                    cmd = 'gunzip *.gz'
                    sp.run(cmd, shell=True)
                else:
                    print('# No output found for {}. Skipping..'.format(name))
                    continue
    print('# Finished')

if __name__ == "__main__":
    main()
