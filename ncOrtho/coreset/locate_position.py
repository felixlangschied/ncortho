
def vprint(s, verbose):
    if verbose:
        print(s, flush=True)


def no_synteny_possible(mirna, chrom, start, end, reference):
    firstgene_info = reference[chrom][1]
    firstgene_start = int(firstgene_info[1])
    lastgene_info = reference[chrom][len(reference[chrom])]
    lastgene_end = int(lastgene_info[2])
    # case 1): there is no protein-coding gene on the same contig as the miRNA,
    # so there can be no neighbors (should only occur in highly fragmented
    # assemblies)
    if chrom not in reference.keys():
        print(
            f'WARNING: No protein-coding genes found on contig "{chrom}". '
            'Make sure that the contig identifiers of the miRNA input file '
            'match the ones in the reference annotation file'
        )
        return True
    # case 2): miRNA is located left of the first gene and hence has no left
    # neighbor, the first gene is therefore by default the right neighbor

    elif end < firstgene_start:
        print(
            f'There is no left neighbor of {mirna}, since it is located at the start of contig {chrom}\n'
            f'{firstgene_info[0]} is the right neighbor of {mirna}'
        )
        return True
    # case 3): miRNA is located right to the last gene, so the last gene is the
    # left neighbor and there cannot be a right neighbor
    elif start > lastgene_end:
        print(
            f'There is no right neighbor of {mirna}, since it is located at the end of contig {chrom}\n'
            f'{lastgene_info[0]} is the left neighbor of {mirna}')
        return True
    else:
        return False


def ortho_search(r_gene, core_dict, v):
    """
    Try to find the ortholog for a given reference gene in a core set species

    Parameters
    ----------
    r_gene:     Gene of Reference Species
    ortho_dict: {'corespecies': {r_gene: [ortho1], ...}, ...}

    Returns:    {'corespecies': [coreortho1], ...}
    -------

    """
    coreorthologs_per_species = {}
    for core_taxon, orthodict in core_dict.items():
        if r_gene in orthodict:
            core_orthologs = orthodict[r_gene]
            coreorthologs_per_species[core_taxon] = core_orthologs
            vprint(f'{" ".join(core_orthologs)} identified as ortholog(s) to {r_gene} in {core_taxon}', v)
        else:
            vprint(f'No ortholog found for {r_gene} in {core_taxon}', v)
    return coreorthologs_per_species


def neighbor_search(leftgene, rightgene, core_dict, gene_position, chromdict, no_next_orthos, v):
    def find_ortho(genename, position, chromd, orthologs, next_orthos, typ, v):
        if genename in orthologs:
            core_orthologs = orthologs[genename]
            vprint(f'{" ".join(core_orthologs)} identified as ortholog(s) to {genename} in {core_taxon}', v)
            return core_orthologs
        else:
            for index in range(1, next_orthos):
                if typ == 'left':
                    nextgeneinfo = chromd[position - index]
                else:
                    nextgeneinfo = chromd[position + index]
                nextgene = nextgeneinfo[0]
                if nextgene in orthologs:
                    core_orthologs = orthologs[nextgene]
                    vprint(f'No  ortholog(s) to {genename} in {core_taxon}', v)
                    vprint(f'{" ".join(core_orthologs)} identified as ortholog(s) of {nextgene}, '
                           f'which is in distance of '
                           f'{index} genes to {genename} in {core_taxon}', v)
                    return core_orthologs

    coreorthologs_per_species = {}
    for core_taxon, orthodict in core_dict.items():
        leftorthos = find_ortho(leftgene, gene_position, chromdict, core_dict[core_taxon], no_next_orthos, 'left', v)
        rightorthos = find_ortho(rightgene, gene_position + 1, chromdict, core_dict[core_taxon], no_next_orthos, 'right', v)
        if leftorthos and rightorthos:
            coreorthologs_per_species[core_taxon] = (leftorthos, rightorthos)
    return coreorthologs_per_species


def categorize_mirna_position(
        mirna, mirna_chrom, mirna_start, mirna_end, mirna_strand, reference, all_orthologs, no_add_orthos, v
):
    """

    Parameters
    ----------
    reference:      Dictionary of reference species annotation as returned by GFF3-Parser
    all_orthologs:  Dictionary of Orthologs as returned by read_pairwise_orthologs
    no_add_orthos:  Number of additional genes up- or downstream that may also serve as syntenic anchors

    Returns

    -------

    """
    syntenytype = None
    solved = False
    core_orthologs = []
    if no_synteny_possible(mirna, mirna_chrom, mirna_start, mirna_end, reference):
        return solved

    # case 4): miRNA is located either between two genes or overlapping with (an
    # intron of) a gene, either on the same or the opposite strand
    ###############################################################################
    for position, geneinfo in reference[mirna_chrom].items():
        gene, gene_start, gene_end, gene_strand = geneinfo
        if mirna_start >= gene_start and mirna_end <= gene_end:
            # case 4.1): miRNA inside gene
            if mirna_strand == gene_strand:
                syntenytype = 'inside'
                solved = True
                vprint(f'{mirna} is located inside the gene {gene}', v)
            # case 4.2): miRNA opposite of gene
            else:
                syntenytype = 'opposite'
                solved = True
                vprint(f'{mirna} is located opposite of the gene {gene}', v)
            core_orthologs = ortho_search(gene, all_orthologs, v)

        # case 4.3): miRNA between genes
        else:
            if position == len(reference[mirna_chrom]):
                continue
            rightneighborinfo = reference[mirna_chrom][position + 1]
            rightneighbor = rightneighborinfo[0]
            rn_start = rightneighborinfo[1]
            if gene_end < mirna_start and rn_start > mirna_end:
                syntenytype = 'in-between'
                solved = True
                vprint(f'{gene} is the left neighbor of {mirna}, {rightneighbor} is the right neighbor', v)
                core_orthologs = neighbor_search(
                    gene, rightneighbor, all_orthologs, position, reference[mirna_chrom], no_add_orthos, v
                )
    if solved:
        return syntenytype, core_orthologs
    else:
        vprint(f'Unable to resolve synteny for {mirna}', v)
        return None, None
