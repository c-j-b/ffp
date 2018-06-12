import sys
from optparse import OptionParser
import re
from scipy import stats, spatial, misc
import subprocess
import os
import time
import numpy as np
from itertools import combinations
import tree_recon


def read_fasta(fasta_fp, seq_num, verbose):
    """
    Used to read in first dictionary and normalize
    :param fasta_fp: fasta input file
    :param seq_num: the number of sequences being compared
    :param verbose: prints  stuff
    :return: The dictionary of the input alignments, where key is the kmer sequences and the value is an array of the freq counts
    """
    fasta_lines = open(fasta_fp).readlines()
    sequences_dct = dict()
    total = 0
    if re.search(">", fasta_lines[0]) is None:
        print("This is not a valid fasta format")
        sys.exit(1)
    for line in fasta_lines:
        if line.strip():
            if re.search(">", line):
                cnt = line.replace(">", "").replace("\n", "").strip()
                total += int(cnt)
            else:
                seq = line.strip("\n")
                cnt_array = [0] * seq_num  # create empty array
                cnt_array[0] = int(cnt)
                sequences_dct[seq] = cnt_array
    # normalize frequency counts
    if verbose:
        print('Normalizing sequence 0')
    for key, value in sequences_dct.items():
        tmp_array = sequences_dct[key]
        tmp_array[0] = value[0]/total
        sequences_dct[key] = tmp_array
    return sequences_dct


def add_fasta(seq_dict, fasta_fp, seq_num, this_seq_num, verbose):
    """
    Used to add second sequences distance value to main dictionary, also normalizes
    :param seq_dict: main dictionary of kmers and frequency counts
    :param this_seq_num: index of this sequence
    :param fasta_fp: second fasta input file
    :param seq_num: total number of sequences being used
    :param verbose: prints stuff
    :return: The dictionary of the input alignments, where key is the kmer sequences and the value is an array of the freq counts
    """
    fasta_lines = open(fasta_fp).readlines()
    total = 0
    if re.search(">", fasta_lines[0]) is None:
        print("This is not a valid fasta format")
        sys.exit(1)
    for line in fasta_lines:
        if line.strip():
            if re.search(">", line):
                cnt = line.replace(">", "").replace("\n", "").strip()
                total += int(cnt)
            else:
                seq = line.strip("\n")
                if seq in seq_dict:
                    # add to array element existing kmer entry
                    tmp_array = seq_dict[seq]
                    tmp_array[this_seq_num] = int(cnt)
                    seq_dict[seq] = tmp_array
                else:
                    # add new kmer
                    cnt_array = [0] * seq_num  # create empty array of max length
                    cnt_array[this_seq_num] = int(cnt)
                    seq_dict[seq] = cnt_array
    # normalize frequency counts
    if verbose:
        print('Normalizing sequence ' + str(this_seq_num))
    for key, value in seq_dict.items():
        tmp_array = seq_dict[key]
        tmp_array[this_seq_num] = value[this_seq_num]/total
        seq_dict[key] = tmp_array
    return seq_dict


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
        #  Count the k-mers
        if verbose:
            print('Jellyfish> Counting k-mers of ' + str(seq_1_file) + ' and ' + str(seq_2_file)
                  + ' with k-mer length ' + str(options.mer_length) + '...')
        # Count k-mers with jellyfish, then dump counts..twice...this takes awhile
        subprocess.run('jellyfish count -m ' + str(options.mer_length) + ' -s 100M -t 10 -C --bf-size 1G '
                       + str(seq_1_file) + ' -o mer_counts_' + str(options.mer_length) + '.jf',
                       cwd=str(seq_1_base), shell=True, check=True)
        if verbose:
            print('Jellyfish> Done counting ' + str(seq_1_file) + '. Starting dump...')
        subprocess.run('jellyfish dump mer_counts_' + str(options.mer_length) +
                       '.jf > mer_counts_dump_' + str(options.mer_length) + '.fa', cwd=str(seq_1_base),
                       shell=True, check=True)
        if verbose:
            print('Jellyfish> Dump complete.')
    else:
        print('\tSkipped counting of ' + seq_1_file)
    if check_for_mer_counts(seq_2_base) is False or rerun_option is True:
        if verbose:
            print('Jellyfish> Counting k-mers of ' + str(seq_2_file) + '...')
        subprocess.run('jellyfish count -m ' + str(options.mer_length) + ' -s 100M -t 10 -C --bf-size 1G '
                       + str(seq_2_file) + ' -o mer_counts_' + str(options.mer_length) + '.jf',
                       cwd=str(seq_2_base), shell=True, check=True)
        if verbose:
            print('Jellyfish> Done counting ' + str(seq_2_file) + '. Starting dump...')
        subprocess.run('jellyfish dump mer_counts_' + str(options.mer_length) +
                       '.jf > mer_counts_dump_' + str(options.mer_length) + '.fa',
                       cwd=str(seq_2_base), shell=True, check=True)
        if verbose:
            print('Jellyfish> Dump complete.')
    else:
        print('\tSkipped counting of ' + seq_2_file)
    # read dumps into dict
    if verbose:
        print('Reading jellyfish FASTAs into dictionaries...')
    seq_dct = read_fasta(str(seq_1_base + '/mer_counts_dump_' + str(options.mer_length) + '.fa'), 2, verbose)
    if verbose:
        print('Done reading into dictionary 1.')
    seq_dct = add_fasta(seq_dct, str(seq_2_base + '/mer_counts_dump_' + str(options.mer_length) + '.fa'), 2, 1, verbose)
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
    ffp_dist = compute_ffp_dist(seq_dct, _distance_option, verbose)
    # print distance
    print('\t' + str(os.path.splitext(os.path.basename(seq_1_file))[0]) + ' + '
          + str(os.path.splitext(os.path.basename(seq_2_file))[0]) + '\n'
          + '\t' + _distance_option + ' dist: ' + str(ffp_dist))
    # write to output file when given path
    if options.output_fp is not None:
        # Write to file
        with open(output_fp, 'a+') as f:
            f.write(str(os.path.splitext(os.path.basename(seq_1_file))[0]) + ' + '
                    + str(os.path.splitext(os.path.basename(seq_2_file))[0]) + '\n' + '\t'
                    + str(ffp_dist) + '\n')
    return ffp_dist


def compute_ffp_dist(sequences_dct, dist_option, verbose):
    """
    :param sequences_dct: The dictionary of the input alignments,
        where keys are the k-mer sequences and values are the k-mer counts in array form
    :param dist_option: Method used to calculate the distance between the FFPs
        options - js: Jenson-Shannon Divergence
                - e: Euclidean Distance
    :param verbose: prints tons of stuff
    :return: the distance between the two FFPs
    """
    if dist_option == 'euclidean':
        if verbose:
            print('Calculating Euclidean distance...')
        d = euclidean([item[0] for item in list(sequences_dct.values())],
                      [item[1] for item in list(sequences_dct.values())])
    elif dist_option == 'euclidean-squared':
        if verbose:
            print('Calculating squared Euclidean distance...')
        d = euclidean_squared([item[0] for item in list(sequences_dct.values())],
                              [item[1] for item in list(sequences_dct.values())])
    else:  # default is JSD
        if verbose:
            print('Calculating Jenson-Shannon divergence...')
        d = jsd([item[0] for item in list(sequences_dct.values())],
                [item[1] for item in list(sequences_dct.values())])
    return d


def check_for_mer_counts(dir_path):
    """
    :param dir_path: folder to look in for mer_counts.jf file
    :return: True if file is found
    """
    flag_1 = False
    flag_2 = False
    for file in os.listdir(dir_path):
        if file == 'mer_counts_' + str(options.mer_length) + '.jf':
            flag_1 = True
        if file == 'mer_counts_dump_' + str(options.mer_length) + '.fa':
            flag_2 = True
    flag = flag_1 and flag_2
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
    parser.add_option("-o", "--output", dest="output_fp", default=None,
                      help="option, path to the output file", metavar="FILE")
    parser.add_option("-d", "--distance-method", dest="distOption", default='js',
                      help="select distance method: Jenson-Shannon (js), Euclidean(e), Euclidean Squared (e2)")
    parser.add_option("-r", "--rerun", dest="rerun", default=True,
                      help="Recounts all k-mers.")
    parser.add_option("-v", "--verbose", dest="verbose", default=False,
                      help="prints lots of stuff")

    (options, args) = parser.parse_args()
    start_time = time.time()

    if options.rerun in ['true', '1', 't', 'y', 'yes', 'True', 'Yes']:
        rerun_bool = True
    else:
        rerun_bool = False
    if options.verbose in ['true', '1', 't', 'y', 'yes', 'True', 'Yes']:
        verbose_bool = True
    else:
        verbose_bool = False
    if options.output_fp is not None:
        # Write to file
        with open(options.output_fp, 'a+') as f:
            f.write('Distance method: ' + str(options.distOption) + '\n')
    if options.folder_fp is None:
        # just use the two sequence paths
        assert options.seq_fp_1 and options.seq_fp_2 is not None, 'Need sequence paths'
        dist = compute_ffp_with_options(options.seq_fp_1, options.seq_fp_2, rerun_bool, options.output_fp,
                                        options.distOption, verbose_bool)
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
            dist = compute_ffp_with_options(file_1, file_2, rerun_bool, options.output_fp,
                                            options.distOption, verbose_bool)
            try:
                result_dict[species_1][species_2] = dist
            except KeyError:
                result_dict[species_1] = {species_2: dist}
            n += 1
            print("\t--- %s seconds ---" % (time.time() - loop_start_time))
        # Reconstruct the tree
        # mirror dictionary for tree reconstruction
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
    if options.output_fp is not None:
        # Write to file
        with open(options.output_fp, 'a+') as f:
            f.write('Running time: ' + str((time.time() - start_time)/60) + '\n')
