import os
import pyfaidx
try:
    from ncOrtho.coreset.coreset_utils import parse_annotation
except ModuleNotFoundError:
    from ncOrtho.coreset.coreset_utils import parse_annotation


def vprint(s, verbose):
    if verbose:
        print(s, flush=True)


def read_genome(path, coretax, outdir):
    core_gen_dir = f'{outdir}/core_genomes'
    if not os.path.isdir(core_gen_dir):
        os.mkdir(core_gen_dir)
    slink = f'{core_gen_dir}/{coretax}'
    try:
        os.symlink(path, slink)
    except FileExistsError:
        pass
    g = pyfaidx.Fasta(slink)
    return g


def synteny_check(left: list, right: list, orthodict: dict, mgi: int, v):
    seqcol = []
    for left_ortho in left:
        if left_ortho not in orthodict:
            continue
        left_chromosome, left_position = orthodict[left_ortho]
        for right_ortho in right:
            if right_ortho not in orthodict:
                continue
            right_chromosome, right_position = orthodict[right_ortho]
            distance = abs(left_position - right_position)
            if left_chromosome == right_chromosome and distance <= mgi:
                leftstart, leftend, leftstrand = orthodict[left_chromosome][left_position][1:4]
                rightstart, rightend, rightstrand = orthodict[left_chromosome][right_position][1:4]
                #vprint('Synteny fulfilled', v)
                if left_position < right_position:
                    seqcol.append((left_chromosome, leftend, rightstart))
                else:
                    seqcol.append((left_chromosome, rightend, leftstart))
    return seqcol


def analyze_synteny(core_d, mirna_pos, out, idtype, mgi, v):
    """

    Parameters
    ----------
    core_d
    mirna_pos
    out
    idtype
    mgi
    v

    Returns:    {mirid: ['>coretaxon1\n', 'seq1\n', 'coretaxon2\n', 'seq2\n']}
    -------

    """
    synteny_region_collector = {}
    for taxon in core_d:
        print(f'# {taxon}', flush=True)
        vprint(f'# Parsing annotation file for {taxon}', v)

        core_anno_dict = parse_annotation(core_d[taxon]['annotation'], idtype)

        vprint('# Loading genome file', v)
        fasta_path = core_d[taxon]['genome']
        genome = read_genome(fasta_path, taxon, out)
        vprint('# Done', v)

        for mirid, positiontuple in mirna_pos.items():
            synteny_fulfilled = False
            if mirid not in synteny_region_collector:
                synteny_region_collector[mirid] = []
            style, core_ortholog_collection_all_taxa = positiontuple
            if taxon not in core_ortholog_collection_all_taxa:
                vprint(f'No core orthologs found for {mirid} in {taxon}', v)
                continue
            core_ortholog_collection = core_ortholog_collection_all_taxa[taxon]
            if style in ['inside', 'opposite']:
                for ortholog in core_ortholog_collection:
                    if ortholog in core_anno_dict:
                        ortholog_chromosome, ortholog_position = core_anno_dict[ortholog]
                        orthostart, orthoend, orthostrand = core_anno_dict[ortholog_chromosome][ortholog_position][1:]
                        seq = genome[ortholog_chromosome][orthostart-1:orthoend].seq

                        synteny_fulfilled = True
                        synteny_region_collector[mirid].append(f'>{taxon}_0\n')
                        synteny_region_collector[mirid].append(f'{seq}\n')
                        # ToDo: check if opposite needs to be reverse
                        # if style == 'inside':
                        #     seq = genome[ortholog_chromosome][orthostart-orthoend].seq
                        # else:
                        #     seq = genome[ortholog_chromosome][orthostart - orthoend].reverse.complement.seq
            else:
                leftorthos, rightorthos = core_ortholog_collection
                synregs = synteny_check(leftorthos, rightorthos, core_anno_dict, mgi, v)
                if not synregs:
                    continue
                for count, syntenyregion in enumerate(synregs):
                    synchrom, synstart, synend = syntenyregion
                    seq = genome[synchrom][synstart-1:synend].seq

                    synteny_fulfilled = True
                    synteny_region_collector[mirid].append(f'>{taxon}_{count}\n')
                    synteny_region_collector[mirid].append(f'{seq}\n')
            if synteny_fulfilled:
                vprint(f'Synteny fulfilled for {mirid}', v)
            else:
                vprint(f'No syntenic region found for {mirid}', v)
    return synteny_region_collector
