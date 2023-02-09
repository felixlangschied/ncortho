import argparse
import sys
import os
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

import ncOrtho.coreset.coreset_utils as cu


def neighbors(prot, anno, anchordist):
    if prot not in anno:  # pseudogenes are not in parsed annotation file but can occur in ortholog file
        return None
    chrom, position = anno[prot]
    anchors = []
    for index in range(position - anchordist, position + anchordist + 1):
        if index not in anno[chrom]:  # if gene positioned at end of contig, no more anchors can be considered
            continue
        anchor, a_start, a_end, a_strand = anno[chrom][index]
        anchors.append(anchor)
    anchors.remove(prot)
    return anchors


def synteny_distance(block, ref_prot, coreanno, orthodict, anchordist, corespecies):
    """
    returns: Number of conserved anchors per number of anchors
    """
    o = []

    # find position of reference protein orthologs in core species
    reforthos = orthodict[ref_prot]
    for refortho in reforthos:
        if refortho not in coreanno:  # pseudogenes are not in parsed annotation file but can occur in ortholog file
            continue
        conserved_anchors = 0
        chrom, position = coreanno[refortho]

        # find position of anchor protein orthologs in core species
        for anchor in block:

            # find orthologs of anchor proteins
            if anchor not in orthodict:  # synteny not fulfilled if no ortholog in core species
                # c.append(0)
                continue
            anchor_ortho_list = orthodict[anchor]
            for anchor_ortho in anchor_ortho_list:
                if anchor_ortho not in coreanno:  # e.g. for pseudogenes
                    # c.append(0)
                    continue

                a_chrom, a_position = coreanno[anchor_ortho]
                # if distance to reference protein is similar in core species than in reference species, accept
                if a_chrom == chrom and abs(int(position) - int(a_position)) <= anchordist:
                    conserved_anchors += 1
                    break

                # c.append(0)
        # o.append(sum(c) / len(c))
        o.append([corespecies, refortho, conserved_anchors])
    return o


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
        '--max_anchor_dist', metavar='int', type=int,
        help='Number of additional genes to the left and right '
             'of the reference miRNA that are to be considered as syntenic anchors. (Default: 3)',
        nargs='?', const=3, default=3
    )
    optional.add_argument(
        '--idtype', metavar='str', type=str,
        help='Choose the ID in the reference gff file that is '
             'compared to the IDs in the pairwise orthologs file (default:GeneID) '
             'Options: ID, Name, GeneID, gene_id, CDS',
        nargs='?', const='ID=', default='ID'
    )
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

    # parameters
    add_pos_orthos = args.max_anchor_dist
    idtype = args.idtype

    core_dict, ref_paths, all_paths = cu.parse_yaml(args.parameters)
    # check if files exist
    # for cp in all_paths:
    #    if not os.path.isfile(cp):
    #        raise ValueError(f'{cp} does not exist')

    ref_anno_dict = cu.gff_parser(ref_paths['annotation'], 'ID')
    print('# Reading pairwise orthologs', flush=True)
    if idtype == 'CDS':
        ortho_dict = cu.pairwise_orthologs_from_cds(core_dict, ref_paths['annotation'])
    else:
        ortho_dict = cu.read_pairwise_orthologs(core_dict)

    syn_col = []
    ortho_col = []
    print('# Checking synteny conservation', flush=True)
    for corespec, pairwiseorthos in ortho_dict.items():

        core_col = []
        print(corespec)
        ortho_col.append([corespec, len(pairwiseorthos.keys())])

        core_anno_dict = cu.gff_parser(core_dict[corespec]['annotation'], 'ID')

        for refprot in pairwiseorthos.keys():
            syn_block = neighbors(refprot, ref_anno_dict, add_pos_orthos)
            if not syn_block:  # in rare case that reference protein from orthology file is not found in annotation
                continue
            syn_fullfilled_list = synteny_distance(syn_block, refprot, core_anno_dict, pairwiseorthos, add_pos_orthos,
                                                   corespec)
            syn_col.extend(syn_fullfilled_list)

    orth = pd.DataFrame(ortho_col, columns=['species', 'Number orthologs'])
    orth = orth.sort_values(by='Number orthologs', ascending=False)
    o = orth.species.values
    sns.set(rc={'figure.figsize': (6, 6), 'ytick.left': True, 'xtick.bottom': True}, font_scale=1.2, style='whitegrid')
    sns.barplot(data=orth, x='species', y='Number orthologs')
    plt.xticks(rotation=45, ha='right')
    plt.xlabel('')
    plt.tight_layout()
    plt.savefig(f'{args.output}/ortholog_distribution.{args.outformat}')


    plt.clf()
    plt.figure(figsize=(12, 4))
    df = pd.DataFrame(syn_col, columns=['species', 'reference_protein', 'num_conserved_anchors'])
    sns.set(rc={'figure.figsize': (12, 4), 'ytick.left': True, 'xtick.bottom': True}, font_scale=1.2, style='whitegrid')
    ax = sns.histplot(data=df, x='num_conserved_anchors', hue='species', discrete=True, multiple='dodge', hue_order=o,
                      shrink=0.8)
    sns.move_legend(ax, "upper left", bbox_to_anchor=(1, 1))
    plt.xlabel(f'Number of conserved anchors in neighborhood k={add_pos_orthos}')
    plt.tight_layout()
    plt.savefig(f'{args.output}/anchor_conservation.{args.outformat}')
    print(f'# Finished. Output created at: {args.output}')


if __name__ == '__main__':
    main()
