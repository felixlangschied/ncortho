# Construct and calibrate a covariance model for a given ncRNA core set
# alignment in Stockholm format.
# The calibration step is computationally very expensive and should be
# performed on multiple CPU cores.
# TODO: Include license, author etc.

import argparse
import multiprocessing as mp
import os
import subprocess as sp
import sys


class CmConstructor(object):
    
    def __init__(self, alignment, outpath, name, cpu):
        self.alignment = alignment
        self.outpath = outpath
        self.name = name
        self.cpu = cpu
        self.model = '{0}/{1}.cm'.format(outpath, name)
    
    def construct(self):
        print('# Constructing covariance model for {}.'.format(self.name))
        construct_command = (
            'cmbuild -n {0} {1}/{0}.cm {2}'
            .format(self.name, self.outpath, self.alignment)
        )
        sp.call(construct_command, shell=True)
        print('# Finished covariance model construction.')
    
    def calibrate(self):
        print('# Calibrating covariance model for {}.'.format(self.name))
        calibrate_command = (
            'cmcalibrate --cpu {0} {1}'.format(self.cpu, self.model)
        )
        sp.call(calibrate_command, shell=True)
        print('# Finished covariance model calibration.')


def main():

    # Parse command-line arguments.
    parser = argparse.ArgumentParser(
        prog='python createcm.py', description='covariance model construction'
    )
    # "cpu", use maximum number of available CPUs unless specified otherwise.
    parser.add_argument(
        '-c', '--cpu', metavar='int', type=int,
        help='number of CPU cores to use', nargs='?',
        const=mp.cpu_count(), default=mp.cpu_count()
    )
    # Path to the desired output folder.
    parser.add_argument(
        '-o', '--output', metavar='<path>', type=str,
        help='path to the output folder'
    )
    # Path to the input alignment
    parser.add_argument(
        '-a', '--alignment', metavar='<.sto>', type=str,
        help='path to input alignment'
    )

    # Show help when no arguments are added.
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)
    else:
        args = parser.parse_args()
    
    # Check if computer provides the desired number of cores.
    available_cpu = mp.cpu_count()
    if args.cpu > available_cpu:
        raise ValueError('The provided number of CPU cores is higher than the number available on this system')
    else:
        cpu = args.cpu

    # Check if alignment file exists.
    if not os.path.isfile(args.alignment):
        raise ValueError('Invalid path to alignment file')
    else:
        alignment = args.alignment

    output = args.output
    name = alignment.split('/')[-1].split('.')[0]

    # Check if the output folder exists.
    if not os.path.isdir(output):
        mkdir_cmd = 'mkdir -p {}'.format(output)
        sp.call(mkdir_cmd, shell=True)

    # Initiate covariance model construction and calibration.
    cmc = CmConstructor(alignment, output, name, cpu)
    # Construct the model.
    cmc.construct()
    # Calibrate the model.
    cmc.calibrate()


if __name__ == '__main__':
    main()
