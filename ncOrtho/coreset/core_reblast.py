import subprocess as sp
import os
import shutil

from ncOrtho.utils import check_blastdb
from ncOrtho.utils import make_blastndb


def vprint(s, verbose):
    if verbose:
        print(s, flush=True)


def make_alignment(out, mirna, cpu, core):
    alignment = '{}/{}.aln'.format(out, mirna)
    stockholm = '{}/{}.sto'.format(out, mirna)

    # Call T-Coffee for the sequence alignment.
    # print('Building T-Coffee alignment.')
    tc_cmd_1 = (
        f't_coffee -quiet -multi_core={cpu} -special_mode=rcoffee -in {core} '
        f'-output=clustalw_aln -outfile={alignment}'
    )
    sp.call(tc_cmd_1, shell=True)

    # Extend the sequence-based alignment by structural information.
    # Create Stockholm alignment.
    # print('Adding secondary structure to Stockholm format.')
    tc_cmd_2 = (
        f't_coffee -other_pg seq_reformat -in {alignment} -action +add_alifold -output stockholm_aln -out {stockholm}'
    )
    sp.call(tc_cmd_2, shell=True)
    return stockholm


def maximum_blast_bitscore(mirna, seq, blastdb, c, dust):
    # The miRNA precursors can show a low level of complexity, hence it might be
    # required to deactivate the dust filter for the BLAST search.
    bit_check = (
        f'blastn -num_threads {c} -dust {dust} -task megablast -db {blastdb} -outfmt \"6 bitscore\"'
    )
    # print('# Determining reference bit score..')
    ref_bit_cmd = sp.Popen(
        bit_check, shell=True, stdin=sp.PIPE, stdout=sp.PIPE, stderr=sp.STDOUT, encoding='utf8'
    )
    ref_results, err = ref_bit_cmd.communicate(seq)
    if not ref_results:
        raise ValueError(f'WARNING: Reference sequence of {mirna} not found in reference Genome. '
              'Consider turning the dust filter off')
    try:
        ref_bit_score = float(ref_results.split('\n')[0].split('\t')[0])
    except ValueError:  # BLASTn errors are in output not in error
        raise ValueError(ref_results)
    return ref_bit_score


# Perform reciprocal BLAST search and construct Stockholm alignment
def blastsearch(mirna, r_path, o_path, c, dust, v):
    """

    Parameters
    ----------
    m_path  :   Path to file with miRNA information
    r_path  :   Path to reference genome.
    o_path  :   Output Name
    c       :   Number of Threads for BLAST search.
    dust    :   Use dust filter for BLAST searches (yes/no)

    Returns
    -------
    Path to Alignment of core pre-miRNAs in Stockholm format

    """

    mirid, rawchrom, mstart, mend, mstrand, rawseq = mirna[:6]
    mchr = rawchrom.replace('chr', '')
    preseq = rawseq.replace('U', 'T')

    miroutdir = f'{o_path}/{mirid}'
    if not os.path.isdir(miroutdir):
        os.mkdir(miroutdir)
    synteny_regs = f'{miroutdir}/synteny_regions_{mirid}.fa'  # this fasta file is created by the the main() script
    os.chdir(miroutdir)

    print(f'# {mirid}')
    if not os.path.isfile(synteny_regs):
        print(f'Warning: No synteny regions found for {mirid}. Training with reference miRNA only.', flush=True)
        with open(synteny_regs, 'w') as fastah:
            fastah.write(f'>{mirid}\n{preseq}\n')
        stock = make_alignment(miroutdir, mirid, c, synteny_regs)
        return stock

    # check if blastdb of reference genome exists
    fname = '.'.join(r_path.split("/")[-1].split('.')[0:-1])
    ref_blastdb = f'{o_path}/refBLASTdb/{fname}'
    if not check_blastdb(ref_blastdb):
        make_blastndb(r_path, ref_blastdb)

    max_bitscore = maximum_blast_bitscore(mirid, preseq, ref_blastdb, c, dust)

    # Find pre-miRNA candidates in syntenic region'
    if not check_blastdb(synteny_regs):
        make_blastndb(synteny_regs, synteny_regs)
    blastn_cmd = (
        'blastn -num_threads {0} -task blastn -dust {1} -db {2} -outfmt \"6 '
        'sseqid evalue bitscore sseq\"'.format(c, dust, synteny_regs)
    )
    blastn = sp.Popen(
        blastn_cmd, shell=True, stdin=sp.PIPE, stdout=sp.PIPE, stderr=sp.PIPE, encoding='utf8'
    )
    ortholog_candidates, err = blastn.communicate(preseq)
    if err:
        raise err

    # Re-BLAST
    outputcol = {}
    seq_check = set()
    ortholog_candidates_list = ortholog_candidates.split('\n')
    vprint(f'Number of ortholog candidates: {len(ortholog_candidates_list)}', v)
    for hit in ortholog_candidates_list:
        if not hit:
            continue
        core_region, eval, bitscore, sseq = hit.strip().split()
        if float(bitscore) <= max_bitscore * 0.5:
            continue

        degap_seq = sseq.replace('-', '')  # Eliminate gaps in BLAST output
        reblastn_cmd = (
            f'blastn -num_threads {c} -task blastn -dust {dust} -db {ref_blastdb} -outfmt \"6'
            ' sseqid sstart send\"'
        )
        reblastn = sp.Popen(
            reblastn_cmd, shell=True, stdin=sp.PIPE,
            stdout=sp.PIPE, stderr=sp.STDOUT, encoding='utf8'
        )
        reresults, reerr = reblastn.communicate(degap_seq)
        for reblast_res in reresults.split('\n'):
            if not reblast_res:
                continue
            refchrom, refstart, refend = reblast_res.split()

            if refchrom != mchr:
                vprint('Hit on chromosome {}, expected {}'.format(refchrom, mchr), v)
                continue
            if (refstart <= mstart <= refend) or (refstart <= mend <= refend):  # first within second
                vprint('Reciprocity fulfilled.', v)
            elif (mstart <= refstart <= mend) or (mstart <= refend <= mend):  # second within first
                vprint('Reciprocity fulfilled.', v)
            else:
                vprint('Reciprocity unfulfilled.', v)
                continue
            if core_region not in seq_check:  # only train with best hit per syntenic region
                outputcol[core_region] = degap_seq
                seq_check.add(core_region)

    corefile = f'{miroutdir}/core_orthologs_{mirid}.fa'
    with open(corefile, 'w') as outfile:
        outfile.write(f'>reference\n{preseq}\n')
        if outputcol:
            for synteny_region, sequence in outputcol.items():
                outfile.write(f'>{synteny_region}\n{sequence}\n')
        else:
            print(f'Warning: No core orthologs found for {mirid}. Training with reference miRNA only.', flush=True)
    stock = make_alignment(miroutdir, mirid, c, corefile)
    return stock
