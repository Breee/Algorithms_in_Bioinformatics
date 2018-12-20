from pprint import pformat

import numpy
from Bio.SubsMat import MatrixInfo

from logger.log import setup_custom_logger
from utility.utils import Alphabet, check_for_duplicates, \
    parse_directory, parse_fasta_files, split_directories_and_files

LOGGER = setup_custom_logger("nw", logfile="gotoh.log")

import argparse
import os


class Xpgma(object):
    """
    Class which implements the WPGMA/UPGMA algorithm.
    """

    def __init__(self):
        LOGGER.info("Initializing")
        sigma = {"A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y"}
        self.alphabet = Alphabet(sigma)
        # The scoring matrix, which is used to calculate the optimal alignment scores.
        self.cluster_matrix = numpy.zeros(shape=(1, 1))

    def create_distance_matrix(self, alignments):
        pass

    def calculate_guide_tree(self):
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
    if args.substitution_matrix != "PAM" and args.substitution_matrix != "BLOSSUM" and args.substitution_matrix != \
            "NONE":
        LOGGER.critical(
                "UNKNOWN parameter for substitution matrix. Choose BLOSSUM or PAM, falling back to BLOSSUM")
        args.substitution_matrix = MatrixInfo.blosum62
    elif args.substitution_matrix == "BLOSSUM":
        args.substitution_matrix = MatrixInfo.blosum62
    elif args.substitution_matrix == "PAM":
        args.substitution_matrix = MatrixInfo.pam250
    elif args.substitution_matrix == "NONE":
        args.substitution_matrix = None


def run_xpgma():
    sequences = parse_input()
    # perform pairwise sequence alignments

    # init the xpgma
    xpgma = Xpgma()
    # create a distance matrix.

    # calculate the guide tree


def main():
    process_program_arguments()
    run_xpgma()
    exit(1)


if __name__ == '__main__':
    # command line arguments.
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input', type=str, required=True, default='',
                        nargs='+',
                        help='A file, a directory or multiple directories. directories are processed recursively.')
    parser.add_argument('--file-filter', type=str, default='',
                        help='A regex to define a filter, which will be applied to the files parsed.')
    parser.add_argument('-g', '--gap_penalty', type=float, default=11,
                        help='Gap start cost, float value. default is 11')
    parser.add_argument('-ge', '--gap_extension', type=float, default=1,
                        help='Gap extension cost, float value. default is 1')
    parser.add_argument('-m', '--match', type=float, default=1,
                        help='match score, only used if no substitution matrix is given, float value. default is 1')
    parser.add_argument('-mm', '--mismatch', type=float, default=-1,
                        help='match score, only used if no substitution matrix is given, float value. default is -1')
    parser.add_argument('-s', '--substitution_matrix', type=str, default="BLOSSUM",
                        help='Substitution Matrix (BLOSSUM | PAM | NONE) default is BLOSSUM')
    parser.add_argument('-d', '--distance', action='store_true', default=False,
                        help='Calculate DISTANCE instead of SIMILARITY')

    args = parser.parse_args()
    main()
else:
    pass
