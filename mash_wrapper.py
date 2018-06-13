#  Chris Brown
#  cjb008@ucsd.edu

from optparse import OptionParser
from scipy import misc
import subprocess
import os
import time
from itertools import combinations


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


if "__main__" == __name__:
    parser = OptionParser()
    parser.add_option("-f", "--folder", dest="folder_fp",
                      help="path to the folder containing leaf sequences", metavar="FILE")
    parser.add_option("-o", "--output", dest="output_fp", default=None,
                      help="option, path to the output file", metavar="FILE")
    parser.add_option("-l", "--length", dest="len", default=30,
                      help="option, path to the output file")
    (options, args) = parser.parse_args()
    start_time = time.time()
    assert options.folder_fp is not None, 'Need sequence folder path'
    if options.output_fp is not None:
        # Write to file
        with open(options.output_fp, 'a+') as f:
            f.write('Distance method: mash\n')
        folders = sorted([f for f in os.listdir(options.folder_fp) if not f.startswith('.')])
        total_combinations = misc.comb(len(folders), 2)
        n = 1
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
            subprocess.run('mash dist -k ' + str(options.len) + str(file_1) + ' ' + str(file_2) + '>> mash_out.txt',
                           cwd=str(options.folder_fp), shell=True, check=True)
            # need to wait for previous command to return
            print("\t--- %s seconds ---" % (time.time() - loop_start_time))
            n += 1
        ("--- %s minutes ---" % ((time.time() - start_time)/60))
