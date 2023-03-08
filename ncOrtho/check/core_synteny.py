import argparse
import sys
import os
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import yaml


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
                refdict['geneset'] = entry['geneset']
            elif in_type == 'core':
                species = entry.pop('name')
                p_dict[species] = entry
            paths.extend([entry['annotation'], entry['geneset']])

        return p_dict, refdict, paths


def read_ortho_tsv(path):
    orthod = {}
    allhumanids = set()
    with open(path) as fh:
        for line in fh:
            if line.startswith('#'):
                continue
            human_protid, core_protid = line.strip().split()
            if not human_protid in orthod:
                orthod[human_protid] = []
            if not core_protid in orthod:
                orthod[core_protid] = []
            orthod[human_protid].append(core_protid)
            orthod[core_protid].append(human_protid)
            allhumanids.add(human_protid)
    return orthod, allhumanids


def geneorder_from_gff(path, proteinset):
    """
    Read gene order from GFF file. Will only collect order of genes for which an ortholog was searched

    returns:
        orderpercontig = {contig: [gene1, gene2, gene3], ...}
        contigperid = {gene1: contig, ...}
    """
    orderpercontig = {}
    contigperid = {}
    lastcontig = ''
    with open(path) as fh:
        for line in fh:
            if line.startswith('#'):
                continue
            data = line.strip().split('\t')
            if data[2] != 'CDS':
                continue
            contig = data[0]
            if contig != lastcontig:
                orderpercontig[contig] = []
                lastcontig = contig
            protid = data[-1].split('-')[1].split('.')[0]

            if orderpercontig[contig] and orderpercontig[contig][-1] == protid:  # do not add duplicates
                continue
            if protid in proteinset:
                orderpercontig[contig].append(protid)
                contigperid[protid] = contig
    return orderpercontig, contigperid


def geneset_from_file(path):
    if path.split('.')[-1] in ['fasta', 'fa', 'faa']:
        o = set()
        with open(path) as fh:
            for line in fh:
                if line.startswith('>'):
                    o.add(line.split('.')[0].replace('>', ''))
    else:
        with open(path) as fh:
            o = [line.strip() for line in fh if line]
            o = set(o)
    return o


def neighbors(protid, orderpercontig, contigperid, k):
    contig = contigperid[protid]
    orderlist = orderpercontig[contig]
    position = orderlist.index(protid)
    # print(contig, position)
    leftflank = orderlist[position - k:position]
    rightflank = orderlist[position + 1:position + k + 1]
    # print(leftflank)
    # print(rightflank)

    return leftflank, rightflank


def orthomapper(referencelist, orthodict):
    o = []
    for prot in referencelist:
        if prot not in orthodict:
            continue
        o.extend(orthodict[prot])
    return o


def jaccard_similarity(list1, list2):
    intersection = len(list(set(list1).intersection(list2)))
    union = (len(set(list1)) + len(set(list2))) - intersection
    if union == 0:
        return union
    else:
        return float(intersection) / union


def main():
    # Print header
    print('\n' + '#' * 38, flush=True)
    print('###' + ' ' * 32 + '###', flush=True)
    print('###   ncOrtho - core set synteny   ###', flush=True)
    print('###' + ' ' * 32 + '###', flush=True)
    print('#' * 38 + '\n', flush=True)

    # Parse command-line arguments
    # Define global variables
    parser = argparse.ArgumentParser(
        description=(
            'Build Covariance models of reference miRNAs from core set of orthologs.'
        )
    )
    parser._action_groups.pop()
    required = parser.add_argument_group('Required Arguments')
    optional = parser.add_argument_group('Optional Arguments')
    # input file path
    required.add_argument('-p', '--parameters', metavar='<path>', type=str, required=True,
                          help='Path to the parameters file in yaml format')
    # output folder
    required.add_argument(
        '-o', '--output', metavar='<path>', type=str, required=True,
        help='Path for the output folder'
    )
    optional.add_argument(
        '-k', '--max_anchor_dist', metavar='int', type=int,
        help='Number of additional genes to the left and right '
             'of the reference miRNA that are to be considered as syntenic anchors. (Default: 3)',
        nargs='?', const=3, default=3
    )
    # optional.add_argument(
    #     '--idtype', metavar='str', type=str,
    #     help='Choose the ID in the reference gff file that is '
    #          'compared to the IDs in the pairwise orthologs file (default:GeneID) '
    #          'Options: ID, Name, GeneID, CDS',
    #     nargs='?', const='ID=', default='ID'
    # )
    optional.add_argument(
        '--outformat', metavar='str', type=str,
        help='Choose format for output figures, as accepted by matplotlib package (Default: png)',
        nargs='?', const='png', default='png'
    )
    # Show help when no arguments are added.
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)
    else:
        args = parser.parse_args()


    core_paths, ref_paths, all_paths = parse_yaml(args.parameters)
    ref_protset = geneset_from_file(ref_paths['geneset'])
    ref_orderpercontig, ref_contigperid = geneorder_from_gff(ref_paths['annotation'], ref_protset)
    print('# Reading reference annotation done', flush=True)

    syn_col = []
    ortho_col = []
    for corespecies, pathdict in core_paths.items():
        print(corespecies)
        pairwise_orthos, allrefids = read_ortho_tsv(pathdict['orthologs'])
        ortho_col.append([corespecies, len(allrefids)])

        core_protset = geneset_from_file(pathdict['geneset'])
        core_orderpercontig, core_contigperid = geneorder_from_gff(pathdict['annotation'], core_protset)

        for refprot in allrefids:
            # check if an ortholog was found for the protein in the reference species
            if not refprot in pairwise_orthos:
                syn_col.append([corespecies, refprot, None])
                continue

            # find position and neighbors of reference protein
            if refprot not in ref_contigperid:  # rare cases of protein.rep.fa and genomic.gff not matching up
                continue
            ref_left, ref_right = neighbors(refprot, ref_orderpercontig, ref_contigperid, args.max_anchor_dist)

            # map neighboring genes to orthologs in core species
            ref_ortho_left = orthomapper(ref_left, pairwise_orthos)
            ref_ortho_right = orthomapper(ref_right, pairwise_orthos)

            # find position and neighbors of the (co-)ortholog to the reference protein in the core species
            coreproteins = pairwise_orthos[refprot]
            for coreprot in coreproteins:
                core_left, core_right = neighbors(coreprot, core_orderpercontig, core_contigperid, args.max_anchor_dist)

                reflist = ref_ortho_left + ref_ortho_right
                corelist = core_left + core_right

                ji = jaccard_similarity(reflist, corelist)
                syn_col.append([corespecies, refprot, ji])


    sns.set(rc={'figure.figsize': (6, 6), 'ytick.left': True, 'xtick.bottom': True}, font_scale=1.2, style='whitegrid')
    orth = pd.DataFrame(ortho_col, columns=['species', 'Number orthologs'])
    sns.barplot(data=orth, x='species', y='Number orthologs', edgecolor="0.3", linewidth=1.5)
    plt.xticks(rotation=45, ha='right')
    plt.xlabel('')
    plt.tight_layout()
    plt.savefig(f'{args.output}/ortholog_distribution.{args.outformat}')


    plt.figure(figsize=(6, 6))
    sns.set(rc={'figure.figsize': (6, 6), 'ytick.left': True, 'xtick.bottom': True}, font_scale=1.2, style='whitegrid')
    df = pd.DataFrame(syn_col, columns=['species', 'reference_protein', 'jaccard_index'])
    df.to_csv(f'{args.output}/conserved_synteny.tsv', sep='\t', index=False)
    sns.boxplot(data=df, x='species', y='jaccard_index', linewidth=1.5)
    plt.xlabel('')
    # plt.title(f'Flanking gene window size = {max_anchor_dist}')
    plt.ylabel(f'Jaccard Similarity (k = {args.max_anchor_dist})')
    plt.xticks(rotation=45, ha='right')
    plt.tight_layout()
    plt.savefig(f'{args.output}/anchor_conservation.{args.outformat}')

    print(f'# Finished. Output created at: {args.output}')

if __name__ == '__main__':
    main()
