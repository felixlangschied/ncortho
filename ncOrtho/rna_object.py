import tempfile
import subprocess as sp
import sys
from tqdm import tqdm
import os
import logging


# Central class of microRNA objects
class Rna(object):
    def __init__(self, name, chromosome, start, end, strand, seq, bit):
        # miRNA identifier
        self.name = name
        # chromosome that the miRNA is located on
        self.chromosome = chromosome
        # start position of the pre-miRNA
        self.start = int(start)
        # end position of the pre-miRNA
        self.end = int(end)
        # sense (+) or anti-sense (-) strand
        self.strand = strand
        # nucleotide sequence of the pre-miRNA
        self.seq = seq.replace('U', 'T')
        # reference bit score that miRNA receives by its own
        # covariance model
        self.bit = bit


def find_cm_refbit(mirid, seq, model, cpu):
    """
    Obtain the reference bit score for each miRNA by applying it to its own covariance model.
    """
    logger = logging.getLogger('ncortho')
    logger.setLevel(level=logging.DEBUG)
    with tempfile.NamedTemporaryFile(mode='w+') as fp:
        fp.write('>{0}\n{1}'.format(mirid, seq))
        fp.seek(0)
        cms_command = f'cmsearch --cpu {cpu} -E 0.01 --noali {model} {fp.name}'
        res = sp.run(cms_command, shell=True, capture_output=True)
    if res.returncode != 0:
        raise sp.SubprocessError(res.stderr.decode("utf-8"))

    for line in res.stdout.decode('utf-8').split('\n'):
        # print(line)
        if line.startswith('  '):
            highest_score = float(line.split()[3])
            return highest_score
    # if no results return topscore 0
    logger.info('Self bit score not applicable, setting threshold to 0')
    return 0.0


def find_phmm_refbit(rnaid, seq, model, cpu):
    logger = logging.getLogger('ncortho')
    logger.setLevel(level=logging.DEBUG)

    call = sp.Popen(
        f'nhmmer -E 0.01 --cpu {cpu} --noali {model} "-"',
        shell=True, stdin=sp.PIPE, stdout=sp.PIPE, stderr=sp.PIPE, encoding='utf8'
        )
    res, err = call.communicate(f'>{rnaid}\n{seq}\n')
    if err:
        raise sp.SubprocessError(err)
    reslist = res.split('\n')
    for i, line in enumerate(reslist):
        if 'Evalue' in line:
            resline = reslist[i+2]
            topscore = float(resline.split()[1])
            return topscore
    # if no results return topscore 0
    logger.info('Self bit score not applicable, setting threshold to 0')
    return 0.0


def rna_maker(mirpath, modeldir, phmm, cpu):
    logger = logging.getLogger('ncortho')
    logger.setLevel(level=logging.DEBUG)
    """
    Parses the miRNA data input file and returns a dictionary of Mirna objects.

    mirpath : Path to file with microRNA data.
    modeldir : Path to directory containing models.

    Returns : {'mirid': Mirna, ...}

    """
    with open(mirpath) as mirna_file:
        rna_data = [
            line.strip().split()[:6] for line in mirna_file if not line.startswith('#')
        ]

    mmdict = {}  # will be the return object
    print('# Calculating reference bit scores')
    for rna in tqdm(rna_data, file=sys.stdout):
        rnaid = rna[0]
        seq = rna[5]

        if phmm:
            model = '{0}/{1}.phmm'.format(modeldir, rnaid)
        else:
            model = '{0}/{1}.cm'.format(modeldir, rnaid)
        if not os.path.isfile(model):
            logger.warning(f'No model found for {rnaid}. Accepted are ".cm" and ".phmm"')
            mmdict[rnaid] = ''
            continue

        if phmm:
            top_score = find_phmm_refbit(rnaid, seq, model, cpu)
        else:
            top_score = find_cm_refbit(rnaid, seq, model, cpu)

        rna.append(top_score)

        # Create output
        mmdict[rnaid] = Rna(*rna)

    return mmdict