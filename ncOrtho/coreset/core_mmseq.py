import os
import subprocess as sp
from glob import glob
import sys
from tempfile import NamedTemporaryFile
from ncOrtho.utils import vprint
from time import time
from datetime import timedelta



def check_mmseq_db(path):
    created_files = glob(f'{path}_h*')
    # print(created_files)
    # index_files = glob(f'{path}.idx*')
    # if len(created_files) > 0 and len(index_files) > 0:
    if len(created_files) > 0:
        return True
    else:
        return False


def make_mmseq_db(mirid, seq, dst):
    with NamedTemporaryFile(mode='w+') as fp:
        fp.write(f'>{mirid}\n{seq}\n')
        fp.seek(0)
        cmd = f'mmseqs createdb {fp.name} {dst} --dbtype 2 -v 2'
        sp.run(cmd, shell=True, check=True)

    # cmd = f'mmseqs createindex {src} {dst} --search-type 3 --mask 0 --split 8'
    # sp.run(cmd, shell=True, check=True)



def mmseq_easy_search(mirid, seq, targetfile, outfile, workingdir, c):
    with NamedTemporaryFile(mode='w+') as fp:
        fp.write(f'>{mirid}\n{seq}\n')
        fp.seek(0)
        cmd = (
            f'mmseqs easy-search {fp.name} {targetfile} {outfile} {workingdir} '
            f'--dbtype 2 --search-type 3 --mask 0 -v 3 --db-load-mode 2 --threads {c} '
            f'--format-output target,tstart,tend,evalue,bits,tseq'
        )
        # print(cmd)
        sp.run(cmd, shell=True, check=True, encoding='utf8')

    with open('{miroutdir}/ortholog_candidates.m8') as fh:
        reslines = [line.strip() for line in fh if line]

    return reslines


def mmseq_search(mirid, seq, targetfile, miroutdir, outfile, workingdir, c):
    make_mmseq_db(mirid, seq, f'{miroutdir}/mirid_mmseqdb')
    cmd = (
        f'mmseqs search {miroutdir}/mirid_mmseqdb {targetfile} {outfile} {workingdir} '
        f'--search-type 3 --mask 0 -v 3 --db-load-mode 2 --threads {c} '
        f'--split-mode 0'

    )
    sp.run(cmd, shell=True, check=True, encoding='utf8')


def find_maximum_bitscore(mirid, seq, reference_database, miroutdir, workingdir, c):
    outfile = f'{miroutdir}/maximum_bitscore'
    mmseq_search(
        mirid, seq, reference_database, miroutdir, outfile, workingdir, c
    )
    print(outfile)
    sys.exit()
    if not reference_hits:
        print(
            f'WARNING: Reference sequence of {mirid} not found in reference Genome. Setting maximum bitscore to 0',
            flush=True
        )
        maximum_bitscore = 0.0
    else:
        best_reference_hit = reference_hits[0]
        maximum_bitscore = best_reference_hit[-2]
    return maximum_bitscore


def coreorthologs_from_mmseq(mirna, refdb, r_path, o_path, c, dust, v):
    mirid, mchr, mstart, mend, mstrand, seq = mirna[:6]
    mchr = mchr.replace('chr', '')
    seq = seq.replace('U', 'T').replace('-', '')
    miroutdir = f'{o_path}/{mirid}'

    synteny_regs = f'{miroutdir}/synteny_regions_{mirid}.fa'  # create by main script

    print(f'# {mirid}', flush=True)
    if not os.path.isfile(synteny_regs):
        print(f'Warning: No synteny regions found for {mirid}. Training with reference miRNA only.', flush=True)
        corefile = f'{miroutdir}/core_orthologs_{mirid}.fa'
        with open(corefile, 'w') as fastah:
            fastah.write(f'>{mirid}\n{seq}\n')
        return corefile

    # check if mmseqdb of reference genome exists
    if refdb:
        if not check_mmseq_db(refdb):
            raise ValueError(f"Cannot find reference MMseq2 database at: {refdb}")
        reference_database = refdb
        dir_path = os.path.realpath(os.path.dirname(refdb))
        workingdir = f'{dir_path}/tmp'
    else:
        if not os.path.isdir(f'{o_path}/reference_database'):
            os.mkdir(f'{o_path}/reference_database')
        reference_database = f'{o_path}/reference_database/refgenome'
        workingdir = f'{o_path}/reference_database/tmp'
        if not os.path.isdir(workingdir):
            os.mkdir(workingdir)
        print('Creating MMseq database. Consider using this path with the "--refdb" option in future or parellel runs:', flush=True)
        print(reference_database, flush=True)

        if not check_mmseq_db(refdb):
            cmd = f'mmseqs createdb {r_path} {reference_database} --dbtype 2 -v 2 --compressed 1'
            sp.run(cmd, shell=True, check=True)

            cmd = f'mmseqs createindex {reference_database} {workingdir} --search-type 3 -v 2'
            sp.run(cmd, shell=True, check=True)

    vprint('Finding maximum bitscore', v)
    st = time()
    maximum_bitscore = find_maximum_bitscore(
        mirid, seq, reference_database, miroutdir, workingdir, c
    )
    et = time()
    td = timedelta(seconds=et - st)
    print('Time used for search in reference:', td)
    print(maximum_bitscore)

    vprint('Search for ortholog candidates', v)
    ortholog_candidates = mmseq_easy_search(
        mirid, seq, synteny_regs, miroutdir, f'{miroutdir}/ortholog_candidates.m8', workingdir, c
    )

    outputcol = {}
    seq_check = set()
    vprint(f'Number of ortholog candidates: {len(ortholog_candidates)}', v)


    sys.exit()






