"""
ncOrtho - Targeted ortholog search for miRNAs
Copyright (C) 2021 Felix Langschied

ncOrtho is a free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

ncOrtho is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with ncOrtho.  If not, see <http://www.gnu.org/licenses/>.
"""

# Modules import

# Python
import argparse
import multiprocessing as mp
import os
import subprocess as sp
import sys
import glob

# Internal ncOrtho modules
from blastparser import BlastParser
from blastparser import ReBlastParser
from genparser import GenomeParser
from cmsearch import cmsearcher

###############################################################################


# Central class of microRNA objects
class Mirna(object):
    def __init__(self, name, chromosome, start, end, strand, pre, mature, bit):
        # miRNA identifier
        self.name = name
        # chromosome that the miRNA is located on
        if chromosome.startswith('chr'):
            self.chromosome = chromosome.split('chr')[1]
        else:
            self.chromosome = chromosome
        # start position of the pre-miRNA
        self.start = int(start)
        # end position of the pre-miRNA
        self.end = int(end)
        # sense (+) or anti-sense (-) strand
        self.strand = strand
        # nucleotide sequence of the pre-miRNA
        self.pre = pre.replace('U', 'T')
        # nucleotide sequence of the mature miRNA
        self.mature = mature.replace('U', 'T')
        # reference bit score that miRNA receives by its own
        # covariance model
        self.bit = bit


def mirna_maker(mirpath, cmpath, output, msl):
    """
    Parses the miRNA data input file and returns a dictionary of Mirna objects.

    Parameters
    ----------
    mirpath : STR
        Path to file with microRNA data.
    cmpath : STR
        Path to covariance models.
    output : STR
        Path for writing temporary files.
    msl : FLOAT
        Length filter.

    Returns
    -------
    mmdict : DICT
        Dictionary with Mirna objects.

    """
    
    mmdict = {} # will be the return object
    
    with open(mirpath) as mirna_file:
        mirna_data = [
            line.strip().split() for line in mirna_file
            if not line.startswith('#')
        ]

    for mirna in mirna_data:
        mirid = mirna[0]
        seq = mirna[5]
        query = '{0}/{1}.fa'.format(output, mirid)
        model = '{0}/{1}.cm'.format(cmpath, mirid)
        # Check if the covariance model even exists, otherwise skip to
        # the next miRNA.
        if not os.path.isfile(model):
            print('# No covariance model found for {}'.format(mirid))
            continue
        # Check if the output folder exists, otherwise create it.
        if not os.path.isdir('{}'.format(output)):
            mkdir = 'mkdir -p {}'.format(output)
            sp.call(mkdir, shell=True)


        # Obtain the reference bit score for each miRNA by applying it
        # to its own covariance model.
        print('# Calculating reference bit score for {}.'.format(mirid))
        
        # Create a temporary FASTA file with the miRNA sequence as
        # query for external search tool cmsearch to calculate the
        # reference bit score.
        with open(query, 'w') as tmpfile:
            tmpfile.write('>{0}\n{1}'.format(mirid, seq))
        cms_output = '{0}/ref_cmsearch_{1}_tmp.out'.format(output, mirid)
        cms_log = '{0}/ref_cmsearch_{1}.log'.format(output, mirid)
        if not os.path.isfile(cms_output):
            cms_command = (
                'cmsearch -E 0.01 --noali -o {3} --tblout {0} {1} {2}'
                .format(cms_output, model, query, cms_log)
            )
            sp.call(cms_command, shell=True)
        else:
            print('# Found cmsearch results at: {} using those'.format(cms_output))
        with open(cms_output) as cmsfile:
            hits = [
                line.strip().split() for line in cmsfile
                if not line.startswith('#')
            ]
            if hits:
                top_score = float(hits[0][14])
            # In case of any issues occuring in the calculation of the bit
            # score, no specific threshold can be determined. The value will
            # set to zero, which technically turns off the filter.
            else:
                print(
                    '# Warning: Self bit score not applicable, '
                    'setting threshold to 0.'
                )
                top_score = 0.0

        mirna.append(top_score)

        # Remove temporary files.
        for rmv_file in [cms_output, cms_log, query]:
            sp.call('rm {}'.format(rmv_file), shell=True)

        # Create output.
        mmdict[mirna[0]] = Mirna(*mirna)

    return mmdict


# blast_search: Perform a reverse BLAST search in the reference genome for a
# candidate.
# s: cmsearch result
# r: reference genome
# o: output name
# c: number of threads
# blast_search(temp_fasta, reference, blast_output, cpu)
def blast_search(s, r, o, c):
    """

    Parameters
    ----------
    s : cmsearch result
    r : Reference genome.
    o : Output Name.
    c : Number of Threads for BLAST search.

    Returns
    -------
    None.

    """
    # Check if BLAST database already exists, otherwise create it.
    # Database files are ".nhr", ".nin", ".nsq".
    file_extensions = ['nhr', 'nin', 'nsq']
    for fe in file_extensions:
        # checkpath = '{}{}'.format(r, fe)
        files = glob.glob(f'{r}*{fe}')
        if not files:
        # At least one of the BLAST db files is not existent and has to be
        # created.
            db_command = 'makeblastdb -in {} -dbtype nucl'.format(r)
            sp.call(db_command, shell=True)
            break
    blast_command = (
       'blastn -task blastn -db {0} -query {1} '
       '-out {2} -num_threads {3} -outfmt "6 qseqid sseqid pident '
       'length mismatch gapopen qstart qend sstart send evalue bitscore sseq"'.format(r, s, o, c)
    )
    sp.call(blast_command, shell=True)


# write_output: Write a FASTA file containing the accepted orthologs.
# Arguments:
# a: dictionary of accepted hits
# o: path for output
# cm: dictionary of cmsearch output
def write_output(a, o, cm):
    """
    

    Parameters
    ----------
    a : TYPE
        DESCRIPTION.
    o : TYPE
        DESCRIPTION.

    Returns
    -------
    None.

    """
    with open(o, 'w') as outfile:
        for hit in a:
            tup = cm[hit]
            header = '|'.join(list(tup))
            outfile.write('>{0}\n{1}\n'.format(header, a[hit]))

# Allow boolean argument parsing
def str2bool(v):
    if isinstance(v, bool):
        return v
    if v.lower() in ('yes', 'true', 't', 'y', '1'):
        return True
    elif v.lower() in ('no', 'false', 'f', 'n', '0'):
        return False
    else:
        raise argparse.ArgumentTypeError('Boolean value expected.')



# Main function
def main():

    # Print header
    print('\n'+'#'*57)
    print('###'+' '*51+'###')
    print('###   ncOrtho - ortholog search for non-coding RNAs   ###')
    print('###'+' '*51+'###')
    print('#'*57+'\n')
    
    # Parse command-line arguments
    # Define global variables
    parser = argparse.ArgumentParser(
        prog='python ncortho.py', description='ncRNA orthology prediction tool'
    )
    # cpu, use maximum number of available cpus unless specified otherwise
    parser.add_argument(
        '-c', '--cpu', metavar='int', type=int,
        help='number of cpu cores ncOrtho should use', nargs='?',
        const=mp.cpu_count(), default=mp.cpu_count()
    )
    # covariance models folder
    parser.add_argument(
        '-m', '--models', metavar='<path>', type=str,
        help='path to your covariance models'
    )
    # mirna data
    parser.add_argument(
        '-n', '--ncrna', metavar='<path>', type=str,
        help='path to your reference micrornas'
    )
    # output folder
    parser.add_argument(
        '-o', '--output', metavar='<path>', type=str,
        help='path for the output folder'
    )
    # query genome
    parser.add_argument(
        '-q', '--query', metavar='<.fa>', type=str,
        help='path to query genome'
    )
    # reference genome
    parser.add_argument(
        '-r', '--reference', metavar='<.fa>', type=str,
        help='path to reference genome'
    )
    # bit score cutoff for cmsearch hits
    parser.add_argument(
        '-t', '--cutoff', metavar='float', type=float,
        help='cmsearch bit score cutoff', nargs='?', const=0.5, default=0.5
    )
    # length filter to prevent short hits
    parser.add_argument(
        '-l', '--msl', metavar='float', type=float,
        help='hit length filter', nargs='?', const=0.7, default=0.7
    )
    parser.add_argument(
        '--heuristic', type=str2bool, metavar='True/False', nargs='?', const=True, default=True,
        help=(
            'Perform a BLAST search of the reference miRNA in the query genome to identify '
            'candidate regions for the cmsearch. Majorly improves speed.'
        )
    )
    parser.add_argument(
        '-x', '--cleanup', type=str2bool, metavar='True/False', nargs='?', const=True, default=True,
        help=(
            'Cleanup temporary files.'
        )
    )
    # check Co-ortholog-ref
    parser.add_argument(
        '--checkCoorthologsRef', type=str2bool, metavar='True/False', nargs='?', const=False, default=False,
        help=(
            'If the re-blast does not identify the original reference miRNA sequence as best hit,'
            'ncOrtho will check whether the best blast '
            'hit is likely a co-ortholog of the reference miRNA relative to the search taxon. '
            'NOTE: Setting this flag will substantially increase'
            'the sensitivity of HaMStR but most likely affect also the specificity, '
            'especially when the search taxon is evolutionarily only very'
            'distantly related to the reference taxon'
        )
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
        print(
            '# Error: The provided number of CPU cores is higher than the '
            'number available on this system. Exiting...'
        )
        sys.exit(1)
    else:
        cpu = args.cpu

    ###TODO: include checks for validity of arguments
    #os.getcwd()
    #os.chdir(path)
    #os.path.exists(path)
    #os.path.isfile(path)
    #os.path.isdir(path)

    mirnas = args.ncrna
    models = args.models
    output = args.output
    query = args.query
    reference = args.reference
    cm_cutoff = args.cutoff
    checkCoorthref = args.checkCoorthologsRef
    cleanup = args.cleanup
    heuristic = args.heuristic
    msl = args.msl
    #blast_cutoff = args.blastc


    # Create miRNA objects from the list of input miRNAs.
    mirna_dict = mirna_maker(mirnas, models, output, msl)

    # Print out query
    print('### Starting ncOrtho run for {}'.format(query))

    # Identify ortholog candidates.
    for mir_data in mirna_dict:
        sys.stdout.flush()
        mirna = mirna_dict[mir_data]
        mirna_id = mirna.name
        if heuristic:
            outdir = '{}'.format(output)
        else:
            outdir = '{}/{}'.format(output, mirna_id)
        # Create output folder, if not existent.
        if not os.path.isdir(outdir):
            os.mkdir(outdir)
        # start cmsearch
        cm_results = cmsearcher(mirna, cm_cutoff, cpu, msl, models, query, output,  cleanup, heuristic)

        # print('cm_results:')
        # print(cm_results)
        # Extract sequences for candidate hits (if any were found).
        if not cm_results:
            print('# No hits found for {}.\n'.format(mirna_id))
            continue
        else:
            gp = GenomeParser(query, cm_results.values())
            candidates = gp.extract_sequences()
            nr_candidates = len(candidates)
            if nr_candidates == 1:
                print(
                    '\n# Covariance model search successful, found 1 '
                    'ortholog candidate.\n'
                )
            else:
                print(
                    '\n# Covariance model search successful, found {} '
                    'ortholog candidates.\n'
                    .format(nr_candidates)
                )
            print('# Evaluating candidates.\n')        
        
        # Perform reverse BLAST test to verify candidates, stored in
        # a list (accepted_hits).
        reblast_hits = {}
        for candidate in candidates:
            sequence = candidates[candidate]
            temp_fasta = '{0}/{1}.fa'.format(outdir, candidate)
            with open(temp_fasta, 'w') as tempfile:
                tempfile.write('>{0}\n{1}\n'.format(candidate, sequence))
            blast_output = '{0}/blast_{1}.out'.format(outdir, candidate)
            print('# Starting reverse blast for {}'.format(candidate))
            # blast_search will write results to the temp_fasta file
            blast_search(temp_fasta, reference, blast_output, cpu)
            # BlastParser will read the temp_fasta file
            bp = BlastParser(mirna, blast_output, msl)
            # the parse_blast_output function will return True if the candidate is accepted
            if bp.evaluate_besthit():
                reblast_hits[candidate] = sequence
            elif checkCoorthref:
                print('Best hit differs from reference sequence! Doing further checks\n')
                # coorth_out = '{0}/{1}_coorth'.format(outdir, candidate)
                if bp.check_coortholog_ref(sequence, outdir):
                    reblast_hits[candidate] = sequence
            if cleanup:
                os.remove(temp_fasta)
                os.remove(blast_output)

        # Write output file if at least one candidate got accepted.
        if reblast_hits:
            nr_accepted = len(reblast_hits)
            if nr_accepted == 1:
                print('# ncOrtho found 1 verified ortholog.\n')
                out_dict = reblast_hits
            else:
                print(
                    '# ncOrtho found {} potential co-orthologs.\n'
                    .format(nr_accepted)
                )
                print('# Evaluating distance between candidates to verify co-orthologs')
                rbp = ReBlastParser(mirna, reblast_hits)
                out_dict = rbp.verify_coorthologs(outdir)
                if len(out_dict) == 1:
                    print('ncOrtho found 1 verified ortholog')
                else:
                    print(
                        '# ncOrtho found {} verified co-orthologs.\n'
                            .format(len(out_dict))
                    )

            print('# Writing output of accepted candidates.\n')
            outpath = '{0}/{1}_orthologs.fa'.format(outdir, mirna_id)
            write_output(out_dict, outpath, cm_results)
            print('# Finished writing output.\n')
        else:
            print(
                '# None of the candidates for {} could be verified.\n'
                .format(mirna_id)
            )
        print('# Finished ortholog search for {}.'.format(mirna_id))

if __name__ == "__main__":
    main()
