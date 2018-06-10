import sys
from optparse import OptionParser
import re
from scipy import stats, spatial, misc
import subprocess
import os
import time
import pandas as pd
import numpy as np
from itertools import combinations
import tree_recon


def read_fasta(fasta_fp):
    """
    :param fasta_fp: fasta input file
    :return: The dictionary of the input alignments, where keys are the kmer sequences and values are the kmer counts
    """
    fasta_lines = open(fasta_fp).readlines()
    sequences_dct = dict()
    if re.search(">", fasta_lines[0]) is None:
        print("This is not a valid fasta format")
        sys.exit(1)
    for line in fasta_lines:
        if line.strip():
            if re.search(">", line):
                cnt = line.replace(">", "").replace("\n", "").strip()
            else:
                seq = line.strip("\n")
                sequences_dct[seq] = int(cnt)
    return sequences_dct


def normalize_counts(d):
    """
    :param d: dictionary with keys as sequences and counts as values
    :return: dictionary with normalized counts
    """
    total = sum(d.values())
    for key, value in d.items():
        d[key] = value / total
    return d


def euclidean(p, q):
    """
    :param p: normalized array like structure
    :param q: normalized array like structure
    :return: Euclidean distance between p & q
    """
    p = np.asarray(p)
    q = np.asarray(q)
    return spatial.distance.euclidean(p, q)


def euclidean_squared(p, q):
    """
    :param p: normalized array like structure
    :param q: normalized array like structure
    :return: Squared euclidean distance between p & q
    """
    p = np.asarray(p)
    q = np.asarray(q)
    return spatial.distance.sqeuclidean(p, q)


def jsd(p, q):
    """
    :param p: normalized array like structure
    :param q: normalized array like structure
    :return: Jenson-Shannon divergence
    """
    p = np.asarray(p)
    q = np.asarray(q)
    m = 0.5 * (p + q)
    return 0.5 * (stats.entropy(p, m) + stats.entropy(q, m))


def compute_ffp_with_options(seq_1_fp, seq_2_fp, rerun_option, output_fp, distance_option, verbose):
    # split full path to get file name and path
    seq_1_base, seq_1_file = os.path.split(seq_1_fp)
    seq_2_base, seq_2_file = os.path.split(seq_2_fp)

    # Automatically skip directories with previous counts unless rerun is enabled
    if check_for_mer_counts(seq_1_base) is False or rerun_option is True:
        #  count the k-mers
        if verbose:
            print('Jellyfish> Counting k-mers of ' + str(seq_1_file) + ' and ' + str(seq_2_file)
                  + ' with k-mer length ' + str(options.mer_length) + '...')
        # count k-mers with jellyfish, then dump counts..twice...this takes awhile
        subprocess.run('jellyfish count -m ' + str(options.mer_length) + ' -s 100M -t 10 -C --bf-size 1G '
                       + str(seq_1_file), cwd=str(seq_1_base), shell=True, check=True)
        if verbose:
            print('Jellyfish> Done counting ' + str(seq_1_file) + '. Starting dump...')
        subprocess.run('jellyfish dump mer_counts.jf > mer_counts_dump.fa ', cwd=str(seq_1_base), shell=True,
                       check=True)
        if verbose:
            print('Jellyfish> Dump complete.')
    else:
        print('\tSkipped counting of ' + seq_1_file)
    if check_for_mer_counts(seq_2_base) is False or rerun_option is True:
        if verbose:
            print('Jellyfish> Counting k-mers of ' + str(seq_2_file) + '...')
        subprocess.run('jellyfish count -m ' + str(options.mer_length) + ' -s 100M -t 10 -C --bf-size 1G '
                       + str(seq_2_file), cwd=str(seq_2_base), shell=True, check=True)
        if verbose:
            print('Jellyfish> Done counting ' + str(seq_2_file) + '. Starting dump...')
        subprocess.run('jellyfish dump mer_counts.jf > mer_counts_dump.fa ', cwd=str(seq_2_base), shell=True,
                       check=True)
        if verbose:
            print('Jellyfish> Dump complete.')
    else:
        print('\tSkipped counting of ' + seq_2_file)
    # read dumps into dict
    if verbose:
        print('Reading jellyfish FASTAs into dictionaries...')
    seq_dct_1 = read_fasta(str(seq_1_base + '/mer_counts_dump.fa'))
    if verbose:
        print('Done reading into dictionary 1.')
    seq_dct_2 = read_fasta(str(seq_2_base + '/mer_counts_dump.fa'))
    if verbose:
        print('Done reading into dictionary 2.')
    # compute distance matrix
    if distance_option == 'e':
        # use euclidean
        _distance_option = "euclidean"
    elif distance_option == 'e2':
        # use euclidean squared
        _distance_option = 'euclidean-squared'
    else:
        # use jenson-shannon
        _distance_option = "jenson-shannon"
    # compute FFP_distance
    ffp_dist = compute_ffp_dist(seq_dct_1, seq_dct_2, _distance_option, verbose)
    # print distance
    print('\t' + str(os.path.splitext(os.path.basename(seq_1_file))[0]) + ' + '
          + str(os.path.splitext(os.path.basename(seq_2_file))[0]) + '\n'
          + '\t\t' + _distance_option + ' dist: ' + str(ffp_dist))
    # write to output file when given path
    if options.output_fp is not None:
        # Write to file
        with open(output_fp, 'a+') as f:
            f.write(str(os.path.splitext(os.path.basename(seq_1_file))[0]) + ' + '
                    + str(os.path.splitext(os.path.basename(seq_2_file))[0]) + '\n' + '\t'
                    + _distance_option + ' dist: ' + str(ffp_dist) + '\n')
    return ffp_dist


def compute_ffp_dist(sequences_dct_1, sequences_dct_2, dist_option, verbose):
    """
    :param sequences_dct_1: The dictionary of the input alignments,
        where keys are the k-mer sequences and values are the k-mer counts normalized
    :param sequences_dct_2: The dictionary of the input alignments,
        where keys are the k-mer sequences and values are the k-mer counts normalized
    :param dist_option: Method used to calculate the distance between the FFPs
        options - js: Jenson-Shannon Divergence
                - e: Euclidean Distance
    :param verbose: prints tons of stuff
    :return: the distance between the two FFPs
    """
    if verbose:
        print('Normalizing FFPs')
    # convert dictionaries to lists of tuples
    norm_seqs_1 = normalize_counts(sequences_dct_1)
    norm_seqs_2 = normalize_counts(sequences_dct_2)
    if verbose:
        print('Aligning k-mers in FFP...')
    # create panda series from dict
    ffp_df = pd.DataFrame.from_records([norm_seqs_1, norm_seqs_2]).fillna(0)
    if verbose:
        print('FFP excerpt:')
        print(ffp_df[ffp_df.columns[:5]])  # print first 5 columns
    if dist_option == 'euclidean':
        if verbose:
            print('Calculating Euclidean distance...')
        d = euclidean(ffp_df.loc[0], ffp_df.loc[1])
    elif dist_option == 'euclidean-squared':
        if verbose:
            print('Calculating squared Euclidean distance...')
        d = euclidean_squared(ffp_df.loc[0], ffp_df.loc[1])
    else:  # default is JSD
        if verbose:
            print('Calculating Jenson-Shannon divergence...')
        d = jsd(ffp_df.loc[0], ffp_df.loc[1])
    return d


def check_for_mer_counts(dir_path):
    """
    :param dir_path: folder to look in for mer_counts.jf file
    :return: True if file is found
    """
    flag = False
    for file in os.listdir(dir_path):
        if file == 'mer_counts.jf':
            flag = True
    flag = flag and check_for_mer_counts_dump(dir_path)
    return flag


def check_for_mer_counts_dump(dir_path):
    """
    :param dir_path: folder to look in for mer_counts_dump.fa file
    :return: True if file is found
    """
    flag = False
    for file in os.listdir(dir_path):
        if file == 'mer_counts_dump.fa':
            flag = True
    return flag


def check_for_fasta(species, base_path):
    """
    :param species: name of species
    :param base_path: base path to species directory
    :return: full path to valid fasta file for species
    """
    fasta_list = []
    full_folder_path = os.path.join(base_path, species)
    for file in os.listdir(full_folder_path):
        if file.endswith('.fq'):
            fasta_list.append(file)
    if len(fasta_list) == 1:
        return os.path.join(base_path, species, fasta_list[0])  # only can have one file
    else:
        return ''


def merge(a, b, path=None):
    """merges b into a"""
    if path is None:
        path = []
    for key in b:
        if key in a:
            if isinstance(a[key], dict) and isinstance(b[key], dict):
                merge(a[key], b[key], path + [str(key)])
            elif a[key] == b[key]:
                pass # same leaf value
            else:
                raise Exception('Conflict at %s' % '.'.join(path + [str(key)]))
        else:
            a[key] = b[key]
    return a


if "__main__" == __name__:
    parser = OptionParser()
    parser.add_option("-f", "--folder", dest="folder_fp",
                      help="path to the folder containing leaf sequences", metavar="FILE")
    parser.add_option("-a", "--sequence1", dest="seq_fp_1",
                      help="path to the FASTA file path", metavar="FILE")
    parser.add_option("-b", "--sequence2", dest="seq_fp_2",
                      help="path to the FASTA file path", metavar="FILE")
    parser.add_option("-l", "--length", dest="mer_length", default=10,
                      help="kmer length to use in jellyfish mer counting")
    parser.add_option("-o", "--output", dest="output_fp",
                      help="option, path to the output file", metavar="FILE")
    parser.add_option("-d", "--distance-method", dest="distOption", default='js',
                      help="select distance method: Jenson-Shannon (js), Euclidean(e), Euclidean Squared (e2)")
    parser.add_option("-r", "--rerun", dest="rerun", default=True,
                      help="Recounts all k-mers.  When false, previously calculated counts are used when available")
    parser.add_option("-v", "--verbose", dest="verbose", default=False,
                      help="prints lots of stuff")

    (options, args) = parser.parse_args()
    start_time = time.time()
    if options.folder_fp is None:
        # just use the two sequence paths
        assert options.seq_fp_1 and options.seq_fp_2 is not None, 'Need sequence paths'
        dist = compute_ffp_with_options(options.seq_fp_1, options.seq_fp_2, options.rerun, options.output_fp,
                                        options.distOption, options.verbose)
    else:
        assert options.folder_fp is not None, 'Need sequence folder path'
        folders = sorted([f for f in os.listdir(options.folder_fp) if not f.startswith('.')])

        total_combinations = misc.comb(len(folders), 2)
        n = 1
        # print('Estimating phylogeny from: ')
        # print(*folders, sep='\n')
        # print('\n')
        result_dict = {}  # Has type result_dict[species_1(str)][species_2(str)] = distance(float)
        # iterate through folder of species building distance dict between all
        for species_1, species_2 in combinations(folders, 2):
            loop_start_time = time.time()
            print('(' + str(n) + '/' + str(total_combinations) + ') '
                  + 'Starting ' + str(species_1) + '  -  ' + str(species_2) + ' distance estimation...')
            file_1 = check_for_fasta(species_1, options.folder_fp)  # check and see if we have a fasta file
            file_2 = check_for_fasta(species_2, options.folder_fp)
            if file_1 == '' or file_2 == '':
                print('\tSkipping combo: ' + str(species_1) + '  -  ' + str(species_2))
                continue
            dist = compute_ffp_with_options(file_1, file_2, options.rerun, options.output_fp,
                                            options.distOption, options.verbose)
            try:
                result_dict[species_1][species_2] = dist
            except KeyError:
                result_dict[species_1] = {species_2: dist}
            n += 1
            print("\t--- %s seconds ---" % (time.time() - loop_start_time))
        print("--- %s minutes ---" % ((time.time() - start_time) / 60))
        # Reconstruct the tree
        # Mirror dictionary
        mirror_result_dict = {}
        for k1, v1 in result_dict.items():  # the basic way
            for k2, v2 in v1.items():
                try:
                    mirror_result_dict[k2][k1] = v2
                except KeyError:
                    mirror_result_dict[k2] = {k1: v2}
        result_dict = merge(result_dict, mirror_result_dict)
        # print(result_dict)
        tree_start_time = time.time()
        print('Reconstructing tree...')
        tree = tree_recon.create_tree(result_dict, folders, 4, 0.3)
        print("--- %s seconds ---" % (time.time() - tree_start_time))
        # write to output file when given path
        if options.output_fp is not None:
            # Write to file
            with open(options.output_fp, 'a+') as f:
                f.write(tree.as_string('newick'))
                f.write(tree.as_ascii_plot())
    print("--- %s minutes ---" % ((time.time() - start_time)/60))
