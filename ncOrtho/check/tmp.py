import argparse
import sys
import os
from collections import Counter
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

import ncOrtho.coreset.coreset_utils as cu

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
             'Options: ID, Name, GeneID, gene_id',
        nargs='?', const='ID=', default='ID'
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
    for cp in all_paths:
        if not os.path.isfile(cp):
            raise ValueError(f'{cp} does not exist')

    print('# Reading pairwise orthologs', flush=True)
    ortho_dict = cu.read_pairwise_orthologs(core_dict)


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


    def synteny_distance(block, ref_prot, coreanno, orthodict, anchordist, corespecies, mgi=3):
        """
        returns: Number of conserved anchors per number of anchors
        """
        o = []

        # find position of reference protein orthologs in core species
        reforthos = orthodict[ref_prot]
        for refortho in reforthos:
            conserved_anchors = 0

            if refortho not in coreanno:  # pseudogenes are not in parsed annotation file but can occur in ortholog file
                continue
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
                    if a_chrom == chrom and abs(int(position) - int(a_position)) <= anchordist + mgi:
                        conserved_anchors += 1
                        break

                    # c.append(0)
            # o.append(sum(c) / len(c))
            o.append([corespecies, refortho, conserved_anchors])
        return o


    print('# Evaluating synteny')
    ref_anno_dict = cu.parse_annotation(ref_paths['annotation'], idtype)
    #syn_col = {}
    syn_col = []
    for corespec, pairwiseorthos in ortho_dict.items():
        #core_col = []
        print(corespec)

        core_anno_dict = cu.parse_annotation(core_dict[corespec]['annotation'], idtype)

        for refprot in pairwiseorthos:
            syn_block = neighbors(refprot, ref_anno_dict, add_pos_orthos)
            if not syn_block:  # in rare case that reference protein from orthology file is not found in annotation
                continue
            syn_fullfilled_list = synteny_distance(syn_block, refprot, core_anno_dict, pairwiseorthos, add_pos_orthos, corespec)
            syn_col.extend(syn_fullfilled_list)

        #syn_col[corespec] = dict(Counter(core_col))

    #print(syn_col)
    df = pd.DataFrame(syn_col, columns=['species', 'reference_protein', 'num_conserved_anchors'])
    #display(df)

    o = ['Gorilla_gorilla', 'Macaca_mulatta', 'Pongo_abelii', 'Nomascus_leucogenys', 'Saimiri_boliviensis',
         'Mus_musculus', 'Cavia_porcellus', 'Canis_familiaris', 'Gallus_gallus']
    sns.set(rc={'figure.figsize': (12, 4), 'ytick.left': True, 'xtick.bottom': True}, font_scale=1.2, style='whitegrid')
    ax = sns.histplot(data=df, x='num_conserved_anchors', hue='species', discrete=True, multiple='dodge', hue_order=o,
                      shrink=0.8)
    sns.move_legend(ax, "upper left", bbox_to_anchor=(1, 1))
    plt.show()



if __name__ == '__main__':
    main()