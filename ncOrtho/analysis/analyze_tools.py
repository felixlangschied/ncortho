import glob
import os
import subprocess as sp
import sys
import shutil
import inspect


def create_overview(results, out):
    # create directory if not esistant
    if not os.path.isdir(out):
        sp.run(f'mkdir -p {out}', shell=True)
    # create overview file or load it if already existant
    outpath = f'{out}/overview.txt'
    ortholog_dict = {}
    print('# Searching for results..')
    files = glob.glob('{}/*/*_orthologs.fa'.format(results))

    with open(outpath, 'w') as of:
        for file in files:
            with open(file, 'r') as infile:
                counter = [line for line in infile if line.startswith('>')]
                count = len(counter)
            of.write(f'{file}\t{count}\n')
            ortholog_dict[file] = count
    print('# Done')
    return ortholog_dict


def make_phyloprofile(overview_dict, map_dict, pp_out):
    print('# Starting to create phyloprofile input')
    # create outdir

    res_dict = {}
    for path in overview_dict:
        count = int(overview_dict[path])
        res_name = path.split('/')[-2]
        mirid = path.split('/')[-1].split('_')[0]

        res_string = '{}#ncbi{}'.format(mirid, map_dict[res_name])
        res_dict[res_string] = count

    with open(pp_out, 'w') as of:
        keys = list(res_dict.keys())
        keys.sort()
        # print(res_C)
        of.write('geneID\tncbiID\torthoID\n')
        for res in keys:
            group, taxid = res.split('#')
            count = res_dict[res]
            for i in range(1, count + 1):
                outstring = '{}\t{}\t{}_c{}\n'.format(group, taxid, group, i)
                of.write(outstring)
    print('# Done')


def extract_representative(data, mirnas, result_path, multi_out):
    # create multifasta files for reach miRNA
    # choosing the top hit ot the cmsearch als represantative ortholog
    print('# Starting to extract representative orthologs in miRNA specific multifasta files')
    if not os.path.isdir(multi_out):
        sp.run(f'mkdir -p {multi_out}', shell=True)

    # load miRNA names
    multifasta_paths = []
    for mirna in mirnas:
        hits = [entry for entry in data if f'/{mirna}_' in entry]
        if not hits:
            print('# Something went run during the extraction of the representative ortholog')
            sys.exit()
        out_file_path = '{}/{}.fa'.format(multi_out, mirna)
        multifasta_paths.append(out_file_path)
        with open(out_file_path, 'w') as of:
            for hit in hits:
                file_2_open = hit.replace('./', f'{result_path}/')
                species = hit.split('/')[-2]

                # if species in spec_to_skip:
                #     continue
                with open(file_2_open, 'r') as readfile:
                    for line in readfile:
                        if not line.startswith('>'):
                            seq = line.strip()
                            break
                of.write(f'>{species}\n')
                of.write(f'{seq}\n')

    print('# Processing of multifasta files complete')
    return multifasta_paths


# aligns multifasta files using mafft
def align_seqs(path_list, align_out, method):
    if not os.path.isdir(align_out):
        sp.run(f'mkdir -p {align_out}', shell=True)
    print('# Starting alignment of {} miRNA sequences using {}'.format(len(path_list), method))
    for path in path_list:
        # output for alignment file
        out_name = path.split('/')[-1].replace('.fa', '.aln')
        out_path = f'{align_out}/{out_name}'
        if method == 'muscle':
            muscle_cmd = f'muscle -in {path} -out {out_path}'
            res = sp.run(muscle_cmd, shell=True, capture_output=True)
            if res.returncode != 0:
                print(res.stderr.decode('utf-8'))
        elif method == 'mafft':
            mafft_cmd = f'linsi {path} > {out_path}'
            res = sp.run(mafft_cmd, shell=True, capture_output=True)
            if res.returncode != 0:
                print(res.stderr.decode('utf-8'))
        else:
            print('ERROR: {} method unknown. Choose between: "mafft", "muscle"'.format(method))
    print('# Alignments finished')


def make_supermatrix(out):
    align_out = '{}/alignments'.format(out)
    tree_out = '{}/supermatrix'.format(out)
    # curr_dir = sys.path[0]
    curr_dir = '/'.join(inspect.getfile(inspect.currentframe()).split('/')[:-1])
    if not os.path.isdir(tree_out):
        sp.run(f'mkdir -p {tree_out}', shell=True)
    print('# Creating supermatrix alignment')
    cmd = (
        'perl {}/concat_alignments_dmp.pl -in {} -out supermatrix.aln'
            .format(curr_dir, align_out)
    )
    res = sp.run(cmd, shell=True, capture_output=True)
    if res.returncode != 0:
        print('ERROR:')
        print(res.stderr.decode('utf-8'))
        sys.exit()
    print('# Done')

    print('# De-gapping alignment')
    cmd = (
        'perl {}/degapper.pl -in {}/supermatrix.aln'.format(curr_dir, out)
    )
    sp.run(cmd, shell=True)
    # move files to output dir
    shutil.move(f'{out}/supermatrix.aln', f'{tree_out}/supermatrix.aln')
    shutil.move(f'{out}/subalignment_positions.txt', f'{tree_out}/subalignment_positions.txt')
    shutil.move(f'{out}/supermatrix.aln.proc', f'{tree_out}/supermatrix.aln.proc')
    print('# Done')
    # return f'{tree_out}/supermatrix.aln'
    return f'{tree_out}/supermatrix.aln.proc'