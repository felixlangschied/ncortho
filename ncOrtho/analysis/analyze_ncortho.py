import os
import subprocess as sp
import sys
import argparse

try:
    from ncOrtho.analysis.analyze_tools import create_overview
    from ncOrtho.analysis.analyze_tools import make_phyloprofile
    from ncOrtho.analysis.analyze_tools import extract_representative
    from ncOrtho.analysis.analyze_tools import align_seqs
    from ncOrtho.analysis.analyze_tools import make_supermatrix
except ImportError:
    from analyze_tools import create_overview
    from analyze_tools import make_phyloprofile
    from analyze_tools import extract_representative
    from analyze_tools import align_seqs
    from analyze_tools import make_supermatrix


def main():
    # Print header
    print('\n'+'#'*39)
    print('###'+' '*33+'###')
    print('###   Analysis of ncOrtho results   ###')
    print('###'+' '*33+'###')
    print('#'*39+'\n')

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
    # mirna data
    required.add_argument(
        '-n', '--ncrna', metavar='<path>', type=str, required=True,
        help='Path to Tab separated file with information about the reference miRNAs. '
             'Can be a subset of all miRNAs for which results are available but must not contain '
             'miRNAs without results in the given directory'
    )
    required.add_argument(
        '-o', '--output', metavar='<path>', type=str,
        help='Path to location where output of analysis should be stored'
    )
    required.add_argument(
        '-m', '--mapping', metavar='<path>', type=str,
        help='Mapping file between NCBI taxonomy id and name of query species/results directory name '
             '(e.g 9606|Homo sapiens)'
    )
    ##################################################################################
    optional.add_argument(
        '--method', metavar='str', type=str, nargs='?', const='mafft', default='mafft',
        help=(
            'Alignment tool to use (Default: MAFFT-linsi) (Options: "mafft", "muscle")'
        )
    )
    optional.add_argument(
        '--skip', metavar='str', type=str, nargs='?', const='', default='',
        help=(
            'Comma separated list of species for which analyses should be skipped'
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
        '--analysis', metavar='str', type=str, nargs='?',
        const='ncOrtho', default='ncOrtho',
        help=(
            'Name of the analysis'
        )
    )
    # Show
    # Show help when no arguments are added.
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)
    else:
        args = parser.parse_args()
    # Parse arguments
    result_path = args.results
    mirna_path = args.ncrna
    map_path = args.mapping
    out_path = args.output
    spec_to_skip = args.skip.split(',')
    tool = args.method
    iqtree_cmd = args.iqtree
    tree_name = args.analysis

    # create output dir if not existing
    if not os.path.isdir(out_path):
        os.makedirs(out_path)

    over_file = '{}/overview.txt'.format(out_path)
    # if overview exists read it
    if os.path.isfile(over_file):
        print('# Ortholog result overview file found')
        ortholog_files = {}
        with open(over_file, 'r') as inf:
            for line in inf:
                key, value = line.strip().split()
                ortholog_files[key] = value

    else:
        # create overview from result directory
        ortholog_files = create_overview(result_path, out_path)

    print('# Reading miRNA info')
    # determine miRNAs to analyse
    mirna_ids = []
    with open(mirna_path, 'r') as inf:
        for line in inf:
            mirna_ids.append(line.strip().split()[0])

    # determine species to analyse
    all_specs = [path.split('/')[-2] for path in ortholog_files]
    all_specs = set(all_specs)
    specs_to_analyse = all_specs - set(spec_to_skip)
    # filter overview dictionary
    filtered_orthologs = {}
    mirnas_found = set()
    # print(mirna_ids)
    for path in ortholog_files:
        outmirna = path.split('/')[-1].replace('_orthologs.fa', '')
        if outmirna in mirna_ids and path.split('/')[-2] in specs_to_analyse:
            mirnas_found.add(path.split('/')[-1].replace('_orthologs.fa', ''))
            filtered_orthologs[path] = ortholog_files[path]
    print('# Done')
    # double check
    if not filtered_orthologs:
        print('# Nothing left after filtering. Exiting..')
        sys.exit()

    # make phyloprofile input
    # 106582|Maylandia zebra|GCF_000238955.4|MAYZE|active
    name_2_id = {}
    with open(map_path, 'r') as mapfile:
        for line in mapfile:
            data = line.strip().split('|')
            taxid = data[0]
            name = data[1].replace(' ', '_')
            name_2_id[name] = taxid
    pp_out = '{}/{}.long'.format(out_path, tree_name)
    make_phyloprofile(filtered_orthologs, name_2_id, pp_out)

    # writing fasta files of representatives to directory in out_path
    multi_out = '{}/multifasta'.format(out_path)
    multi_paths = extract_representative(filtered_orthologs, mirnas_found, result_path, multi_out)


    # align representative multifasta files using MAFFT-linsi
    align_out = '{}/alignments'.format(out_path)
    align_seqs(multi_paths, align_out, tool)

    # call Ingo's concat alignment and degap scripts
    superm_path = make_supermatrix(out_path)

    # calculate tree using iqtree
    print('# Starting tree calculation')
    tree_out = f'{out_path}/{tree_name}'
    if not os.path.isdir(tree_out):
        os.mkdir(tree_out)
    tree_cmd = iqtree_cmd.format(superm_path, tree_out)
    res = sp.run(tree_cmd, shell=True, capture_output=True)
    print(res.stdout.decode('utf-8'))
    if res.returncode != 0:
        print(res.stderr.decode('utf-8'))


if __name__ == '__main__':
    main()