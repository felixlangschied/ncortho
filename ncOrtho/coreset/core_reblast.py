import subprocess as sp
import os
import shutil
from time import time
from datetime import timedelta
from tqdm import tqdm

from ncOrtho.utils import check_blastdb
from ncOrtho.utils import make_blastndb


def vprint(s, verbose):
    if verbose:
        print(s, flush=True)


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
        print(f'WARNING: Reference sequence of {mirna} not found in reference Genome. Setting maximum bitscore to 0', flush=True)
        return 0.0
    try:
        ref_bit_score = float(ref_results.split('\n')[0].split('\t')[0])
    except ValueError:  # BLASTn errors are in output not in error
        print(f'WARNING: Reference sequence of {mirna} not found in reference Genome. Setting maximum bitscore to 0', flush=True)
        return 0.0
    return ref_bit_score


# Perform reciprocal BLAST search and construct Stockholm alignment
def blastsearch(mirna, refdb, r_path, o_path, c, dust, reblastmode, v):
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
    preseq = rawseq.replace('U', 'T').replace('-', '')

    miroutdir = f'{o_path}/{mirid}'
    if not os.path.isdir(miroutdir):
        os.mkdir(miroutdir)
    synteny_regs = f'{miroutdir}/synteny_regions_{mirid}.fa'  # this fasta file is created by the the main() script
    os.chdir(miroutdir)

    print(f'# {mirid}', flush=True)
    if not os.path.isfile(synteny_regs):
        print(f'Warning: No synteny regions found for {mirid}. Training with reference miRNA only.', flush=True)
        corefile = f'{miroutdir}/core_orthologs_{mirid}.fa'
        with open(corefile, 'w') as fastah:
            fastah.write(f'>reference\n{preseq}\n')
        return corefile
        # stock = make_alignment(miroutdir, mirid, c, synteny_regs, coffee)
        # return stock

    # check if blastdb of reference genome exists
    if refdb:
        if check_blastdb(refdb):
            ref_blastdb = refdb
        else:
            raise ValueError(f"Cannot find reference BLAST databse at: {refdb}")
    else:
        fname = '.'.join(r_path.split("/")[-1].split('.')[0:-1])
        ref_blastdb = f'{o_path}/reference_database/{fname}'
        if not check_blastdb(ref_blastdb):
            make_blastndb(r_path, ref_blastdb)

    st = time()
    max_bitscore = maximum_blast_bitscore(mirid, preseq, ref_blastdb, c, dust)
    et = time()
    td = timedelta(seconds=et - st)
    print('Time to find maximum bitscore in reference:', td)

    # Find pre-miRNA candidates in syntenic region'
    st = time()
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
        raise sp.SubprocessError(err)
    et = time()
    td = timedelta(seconds=et - st)
    print('Time to find ortholog candidates:', td)

    # Re-BLAST
    outputcol = {}
    seq_check = set()
    ortholog_candidates_list = ortholog_candidates.split('\n')

    # vprint(f'Number of ortholog candidates: {len(ortholog_candidates_list)}', v)
    for hit in tqdm(ortholog_candidates_list):
        if not hit:
            continue
        core_region, eval, bitscore, sseq = hit.strip().split()
        if core_region in seq_check:
            continue

        if float(bitscore) <= max_bitscore * 0.5:
            continue

        vprint(f'{core_region}, length: {len(sseq)} nt', v)

        degap_seq = sseq.replace('-', '')  # Eliminate gaps in BLAST output

        if reblastmode == 'blastn':
            # sensitive mode for miRNAs
            reblastn_cmd = (
                f'blastn -num_threads {c} -task blastn -dust {dust} -db {ref_blastdb} -outfmt \"6'
                ' sseqid sstart send\"'
            )
        else:
            # efficient mode for lincRNAs
            reblastn_cmd = (
                f'blastn -num_threads {c} -task megablast -dust {dust} -db {ref_blastdb} -outfmt \"6'
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
            # print(reblast_res)
            refchrom, refstart, refend = reblast_res.split()

            if refchrom != mchr:
                # vprint('Hit on chromosome {}, expected {}'.format(refchrom, mchr), v)
                continue
            if (refstart <= mstart <= refend) or (refstart <= mend <= refend):  # first within second
                vprint('Reciprocity fulfilled.', v)
            elif (mstart <= refstart <= mend) or (mstart <= refend <= mend):  # second within first
                vprint('Reciprocity fulfilled.', v)
            else:
                # vprint('Reciprocity unfulfilled.', v)
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

    return corefile
