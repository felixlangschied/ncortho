import glob
import os
import tempfile
import subprocess as sp
import sys
import shutil
import inspect
from collections import Counter
import argparse

"""
Species in spec_to_skip will not be part of multiple sequence alignments
"""


def make_supermatrix(out):
    align_out = '{}/alignments'.format(out)
    tree_out = '{}/supermatrix'.format(out)
    # curr_dir = sys.path[0]
    curr_dir = '/'.join(inspect.getfile(inspect.currentframe()).split('/')[:-1])
    if not os.path.isdir(tree_out):
        os.mkdir(tree_out)
    print('# Creating supermatrix alignment')
    cmd = (
        'perl {}/concat_alignments_dmp.pl -in {} -out supermatrix.aln'
            .format(curr_dir, align_out)
    )
    res = sp.run(cmd, shell=True, capture_output=True)
    if res.returncode != 0:
        print('ERROR:')
        print(res.stderr.decode('utf-8'))
        sys.exit()
    # print('# Done')

    print('# De-gapping alignment')
    cmd = (
        'perl {}/degapper.pl -in {}/supermatrix.aln'.format(curr_dir, out)
    )
    sp.run(cmd, shell=True)
    # move files to output dir
    shutil.move(f'{out}/supermatrix.aln', f'{tree_out}/supermatrix.aln')
    shutil.move(f'{out}/subalignment_positions.txt', f'{tree_out}/subalignment_positions.txt')
    shutil.move(f'{out}/supermatrix.aln.proc', f'{tree_out}/supermatrix.aln.proc')
    # print('# Done')
    # return f'{tree_out}/supermatrix.aln'
    return f'{tree_out}/supermatrix.aln.proc'


def arglistcheck(arg):
    if not arg:
        return ''
    if os.path.isfile(arg):
        with open(arg) as fh:
            arglist = [element.strip() for element in fh.read().split('\n') if element]
    else:
        arglist = [element.strip() for element in arg.split(',') if element]
    return arglist


def parse_matseq(path):
    m2s = {}
    with open(path) as fh:
        for line in fh:
            if line.startswith('#'):
                continue
            mature = line.strip().split('\t')[-1].replace('U', 'T')
            mirid = line.strip().split('\t')[1].split('_')[0]
            if len(mature) < 7:
                raise ValueError('No mature sequence found in "ncrna_file"')
            seed = mature[1:8]
            m2s[mirid] = seed
    return m2s

def main():
    # Print header
    print('\n'+'#'*39)
    print('###'+' '*33+'###')
    print('###   Analysis of ncOrtho results   ###')
    print('###'+' '*33+'###')
    print('#'*39+'\n')
    sys.stdout.flush()

    # Parse command-line arguments
    # Define global variables
    parser = argparse.ArgumentParser(
        description='Analzye ncOrtho results. Gives: PhyloProfile input, supermatrix alignment, species tree.'
    )
    parser._action_groups.pop()
    required = parser.add_argument_group('Required Arguments')
    optional = parser.add_argument_group('Optional Arguments')
    required.add_argument(
        '-r', '--results', metavar='<path>', type=str,
        help='Path to ncOrtho output directory that is to be analyzed'
    )
    required.add_argument(
        '-o', '--output', metavar='<path>', type=str,
        help='Path to location where output of analysis should be stored. '
             'Assumes that result files end with "_orthologs.fa"'
    )
    required.add_argument(
        '-m', '--mapping', metavar='<path>', type=str,
        help='Tab seperated mapping file between NCBI taxonomy id and name of the query species '
             'or more generally the name of the directory that contains the ncOrtho results of each species'
             '(e.g "9606\tHomo_sapiens")'
    )
    ##########################################################################################
    optional.add_argument(
        '--skip', metavar='str', type=str, nargs='?', const='', default='',
        help=(
            'Comma separated list or path to newline seperated list of species for which analyses should be skipped'
        )
    )
    optional.add_argument(
        '--include', metavar='str', type=str, nargs='?', const='', default='',
        help=(
            'Comma separated list or path to newline seperated list of species for which analyses should be performed. '
            'Will skip every species not in list.'
        )
    )
    optional.add_argument(
        '--auto_skip', metavar='float', type=float, nargs='?', const=0.0, default=0.0,
        help=(
            'Exclude species for which less than a fraction of the reference '
            'miRNAs were detected from species tree calculation. '
            '(e.g. --auto_skip 0.5 means that only species in which at least 50%% '
            'of the reference miRNAs were found are used for species tree calculation) (Default: Off)'
        )
    )
    optional.add_argument(
        '--iqtree', metavar='str', type=str, nargs='?',
        const='iqtree -s {} -bb 1000 -alrt 1000 -nt AUTO -redo -pre {}/species_tree',
        default='iqtree -s {} -bb 1000 -alrt 1000 -nt AUTO -redo -pre {}/species_tree',
        help=(
            'Call to iqtree. Do not change curly bracket notation '
            '(Default: iqtree -s {} -bb 1000 -alrt 1000 -nt AUTO -redo -pre {}/{})'
        )
    )
    optional.add_argument(
        '--mirnas', metavar='<path>', type=str, nargs='?', const='', default='',
        help=(
            'Comma separated list or path to newline seperated file with miRNAs that should be used for reconstruction'
        )
    )
    # mirna data
    optional.add_argument(
        '--ncrna_file', metavar='<path>', nargs='?', const='', default='',
        help='Path to tab separated file of reference miRNAs information, as used in ncCreate and ncSearch. '
             'Uses mature sequence column to identify seed conservation.'
    )
    ##########################################################################################
    # Parse arguments
    ##########################################################################################
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)
    else:
        args = parser.parse_args()

    res_dict = args.results
    taxmap = args.mapping
    outdir = args.output

    spec_to_skip = arglistcheck(args.skip)
    spec_include = arglistcheck(args.include)

    overlap = set(spec_to_skip).intersection(set(spec_include))
    if overlap:
        raise ValueError(f'"{overlap}" present in species to skip as well as species to include')


    mirlist = arglistcheck(args.mirnas)

    iqtree_cmd = args.iqtree
    auto_skip = args.auto_skip

    if args.ncrna_file:
        mirid2seed = parse_matseq(args.ncrna_file)
    else:
        mirid2seed = {}
    ##########################################################################################
    # Argument checks
    ##########################################################################################
    if not os.path.isdir(outdir):
        os.mkdir(outdir)
    if auto_skip < 0 or auto_skip > 1:
        raise ValueError(f'--auto_skip must be a value between 0 and 1, not {auto_skip}')

    ###############################################################################################
    # Main
    ###############################################################################################
    # find results
    ortholog_files = glob.glob(f'{res_dict}/*/*_orthologs.fa')
    if not ortholog_files:
        ortholog_files = glob.glob(f'{res_dict}/*.fa')
    # read mapping file
    name_2_taxid = {}
    with open(taxmap, 'r') as mp:
        for line in mp:
            taxid, name = line.strip().split('\t')
            name_2_taxid[name] = taxid

    ortho_dict = {}  # ortho_dict = {mirna: {species: seq}}
    # read all results into one dictionary and write phyloprofile input
    pp_in = f'{outdir}/PhyloProfile.long'
    spec_list = []
    with open(pp_in, 'w') as pp:
        pp.write('geneID\tncbiID\torthoID\n')
        for file in ortholog_files:
            with open(file, 'r') as fh:
                for line in fh:
                    if line.startswith('>'):
                        spec, mirna = line.strip().replace('>', '').split('|')[:2]
                    else:
                        seq = line.strip()
                        # write phyloprofile input
                        taxstr = f'ncbi{name_2_taxid[spec]}'
                        group_str = mirna.split('_')[0]
                        if mirid2seed:
                            seed = mirid2seed[group_str]
                            if seed in seq:
                                seedcon = 1
                            else:
                                seedcon = 0
                            pp.write(f'{group_str}\t{taxstr}\t{mirna}\t{seedcon}\n')
                        else:
                            pp.write(f'{group_str}\t{taxstr}\t{mirna}\n')
                        # initialize dictionary
                        if group_str not in ortho_dict:
                            ortho_dict[group_str] = {}
                        # only save representative (first) hit
                        if spec not in ortho_dict[group_str]:
                            ortho_dict[group_str][spec] = seq
                            spec_list.append(spec)
    print('# Finished writing PhyloProfile input')

    # if miRNAfile option is chosen, prune dictionary
    if mirlist:
        tmpdict = {}
        for mirnaid, value in ortho_dict.items():
            if mirnaid in mirlist:
                tmpdict[mirnaid] = value
        del ortho_dict
        ortho_dict = tmpdict

    # count species and skip species with few orthologs before aligning
    if auto_skip:
        spec_counter = Counter(spec_list)
        max_spec_count = spec_counter.most_common(1)[0][1]
        count_cutoff = max_spec_count * auto_skip
        print(f'# Only species with at least {round(count_cutoff)} identified orthologs are used for tree calculation')
        for species, count in dict(spec_counter).items():
            if count < count_cutoff:
                spec_to_skip.append(species)
        print('# Skipping these species:')
        for spsk in spec_to_skip:
            print(spsk)

    # if no species list is given to include, all species are part of specinclude
    if not spec_include:
        spec_include = list(ortho_dict.keys())

    # make alignments
    print('# Starting alignments')
    align_out = f'{outdir}/alignments'
    if not os.path.isdir(align_out):
        os.mkdir(align_out)
    for mirna in ortho_dict:
        with tempfile.NamedTemporaryFile(mode='w+') as fp:
            for spec, seq in ortho_dict[mirna].items():
                if (
                        spec not in spec_to_skip
                        and spec in spec_include
                ):
                    fp.write(f'>{spec}\n{seq}\n')
            fp.seek(0)
            aln_cmd = f'muscle -in {fp.name} -out {align_out}/{mirna}.aln'
            res = sp.run(aln_cmd, shell=True, capture_output=True)
            if res.returncode != 0:
                print(res.stderr.decode('utf8'))
                sys.exit()

    sys.stdout.flush()
    superm_path = make_supermatrix(outdir)
    print('# Starting tree calculation')
    tree_out = f'{outdir}/species_tree'
    if not os.path.isdir(tree_out):
        os.mkdir(tree_out)
    tree_cmd = iqtree_cmd.format(superm_path, tree_out)
    res = sp.run(tree_cmd, shell=True, capture_output=True)
    print(res.stdout.decode('utf-8'))
    if res.returncode != 0:
        print(res.stderr.decode('utf-8'))

    # write tree that can be used as optional phyloprofile input
    treepath = f'{tree_out}/species_tree.contree'
    pptree = f'{outdir}/PhyloProfile_tree.contree'
    with open(treepath, 'r') as th:
        treestr = th.read().strip()
        for name in name_2_taxid:
            ppid = f'ncbi{name_2_taxid[name]}'
            treestr = treestr.replace(name, ppid)
    with open(pptree, 'w') as otfh:
        otfh.write(treestr)


if __name__ == "__main__":
    main()
