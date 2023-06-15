import glob
import subprocess as sp
import argparse


def check_blastdb(db_path):
    file_extensions = ['nhr', 'nin', 'nsq']
    for fe in file_extensions:
        files = glob.glob(f'{db_path}*{fe}')
        if not files:
            # At least one of the BLAST db files is not existent
            return False
    return True


def make_blastndb(inpath, outpath):
    db_command = 'makeblastdb -in {} -out {} -dbtype nucl'.format(inpath, outpath)
    sp.run(db_command, shell=True, capture_output=True)


def vprint(s, verbose):
    if verbose:
        print(s, flush=True)


def write_output(outlist, path):
    with open(path, 'w') as of:
        for line in outlist:
            of.write(line)


def str2bool(v):
    if isinstance(v, bool):
        return v
    if v.lower() in ('yes', 'true', 't', 'y', '1'):
        return True
    elif v.lower() in ('no', 'false', 'f', 'n', '0'):
        return False
    else:
        raise argparse.ArgumentTypeError('Boolean value expected.')