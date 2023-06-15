try:
    import subprocess as sp
    import logging
    from blastparser import BlastParser
except ModuleNotFoundError:
    from ncOrtho.blastparser import BlastParser


def perform_reblast(sequence, refblast, cpu, outdir, candidate, mirna_data, dust, msl, cleanup, checkCoorthref):
    logger = logging.getLogger('ncortho')
    logger.setLevel(level=logging.DEBUG)

    blast_command = (
        f'blastn -task blastn -db {refblast} -num_threads {cpu} -dust {dust} '
        f'-outfmt "6 sseqid sstart send sstrand bitscore"'
    )
    blast_call = sp.Popen(
        blast_command, shell=True, stdin=sp.PIPE,
        stdout=sp.PIPE, stderr=sp.PIPE, encoding='utf8'
    )
    # run BLAST
    res, err = blast_call.communicate(sequence)
    if err:
        logger.warning(err)
        return ''
        # continue
    if not cleanup:
        blast_output = f'{outdir}/reBLAST_{candidate}.out'
        with open(blast_output, 'w') as bh:
            for line in res:
                bh.write(line)
    # parse BLASTn results
    blast_output = [line.split() for line in res.split('\n')]
    if not blast_output[0]:
        logger.info('No re-BLAST hits')
        return ''

    # BlastParser will read the temp_fasta file
    bp = BlastParser(mirna_data, blast_output, msl)
    # the parse_blast_output function will return True if the candidate is accepted
    if bp.evaluate_besthit():
        logger.info('Found best hit')
        return sequence

    elif checkCoorthref:
        logger.info('Best hit differs from reference sequence! Doing further checks\n')
        # coorth_out = '{0}/{1}_coorth'.format(outdir, candidate)
        if bp.check_coortholog_ref(sequence, outdir):
            return sequence
    else:
        logger.info('Best hit does not overlap with miRNA location')
        return ''