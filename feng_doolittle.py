from pprint import pformat

from logger.log import setup_custom_logger
from utility.utils import check_for_duplicates, \
    parse_directory, parse_fasta_files, split_directories_and_files

LOGGER = setup_custom_logger("got", logfile="gotoh.log")

import argparse
import os


class FengDoolittle(object):
    """
    Class which implements the feng doolittle
    1. Calculate the 38#38 pairwise alignment scores, and convert them to distances.
    2. Use an incremental clustering algorithm (Fitch and Margoliash, 1967 [3]) to construct a tree from the distances.
    3. Traverse the nodes in their order of addition to the tree, repeatedly align the child nodes (sequences or alignments).
    Features of this heuristic:
    - Highest scoring pairwise alignment determines the alignment to two groups.
    - "Once a gap, always a gap": replace gaps in alignments by a neutral character.
    """

    def __init__(self, similarity=True, verbose=False):
        pass


def parse_input():
    fasta_files = []
    # split input into files and directories.
    directories, files = split_directories_and_files(input_list=args.input)
    # check if input files are part of the directories to be checked
    # an  check if directories are subdirectories of other directories.
    directories, files = check_for_duplicates(directories=directories, files=files)
    for file in files:
        fasta_files.append(file)
        # process directories and get fastafiles.
    for dir_name in directories:
        directory_content = parse_directory(dir_name, file_filter=args.file_filter)
        for entry in directory_content:
            os.chdir(entry["directory"])
            if entry["files"] != [] and entry["directory"] != '':
                fasta_files.extend(entry["files"])
    LOGGER.info("Collected the following fasta files:\n %s" % pformat(fasta_files))
    sequences = parse_fasta_files(fasta_files)
    LOGGER.info("Parsed the following sequences:\n %s" % pformat(sequences))
    return sequences


def process_program_arguments():
    if not args.input:
        LOGGER.critical(
                "Error, provide input file/files/directory/directories using -i / --input. -h/--help for all "
                "arguments.")
        exit(1)


def run_feng_doolittle():
    sequences = parse_input()
    if len(sequences) < 2:
        LOGGER.warn("We received not enough sequences. Make sure you called the program correctly.")
        exit(1)
    elif len(sequences) >= 2:
        pass


def main():
    process_program_arguments()
    run_feng_doolittle()
    exit(1)


if __name__ == '__main__':
    # command line arguments.
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input', type=str, required=True, default='',
                        nargs='+',
                        help='A file, a directory or multiple directories. directories are processed recursively.')
    parser.add_argument('--file-filter', type=str, default='',
                        help='A regex to define a filter, which will be applied to the files parsed.')
    parser.add_argument('-v', '--verbose', action='store_true', default=False,
                        help='verbose output.')

    args = parser.parse_args()
    main()
else:
    pass
