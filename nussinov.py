"""
MIT License

Copyright (c) 2018 Julian Löffler (Breee@github)

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
"""

import logging

import numpy as np

from logger.log import setup_custom_logger
from utility.utils import Alphabet, parse_input

LOGGER = setup_custom_logger("nov", logfile="nussinov.log")

import argparse


class Nussinov(object):
    """
    Class which implements the nussinov algorithm.
    """

    def __init__(self, min_loop_length=4, verbose=False):
        LOGGER.info("Initialzing nussinov")
        sigma = {"A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y",
                 "X", "U"}
        self.alphabet = Alphabet(sigma)
        self.matrix = None
        self.structure = []
        self.basepairs = {frozenset(["A", "U"]), frozenset(["G", "C"])}
        self.min_loop_length = min_loop_length
        if verbose:
            LOGGER.level = logging.DEBUG

    def is_basepair(self, letter1, letter2):
        return {letter1, letter2} in self.basepairs

    def structure_to_brackets_coords(self, sequence):
        """
        Convert a structure to bracket representation
        :param sequence: input sequence.
        :return:  bracket representation
        """
        assert self.structure, "structure is empty."
        dot_list = ["." for _ in range(len(sequence))]
        coords_list = []
        for s in self.structure:
            coords = s['coords']
            coords_list.append((min(coords), max(coords)))
            dot_list[min(coords)] = "("
            dot_list[max(coords)] = ")"
        return "".join(dot_list), coords_list

    def traceback(self, sequence, i, j, structure):
        if i < j:
            if self.matrix[i][j] == self.matrix[i + 1][j]:
                self.traceback(sequence, i + 1, j, structure)
            elif self.matrix[i][j] == self.matrix[i][j - 1]:
                self.traceback(sequence, i, j - 1, structure)
            elif self.matrix[i][j] == self.matrix[i + 1][j - 1] + 1 if self.is_basepair(sequence[i],
                                                                                        sequence[j]) else 0:
                structure.append({'coords': (i, j), 'letters': (str(sequence[i]), str(sequence[j]))})
                self.traceback(sequence, i + 1, j - 1, structure)
            else:
                for k in range(i + 1, j):
                    if self.matrix[i][j] == self.matrix[i][k] + self.matrix[k + 1][j]:
                        self.traceback(sequence, i, k, structure)
                        self.traceback(sequence, k + 1, j, structure)
                        break
        return structure

    def run(self, sequence):
        """
        Function which runs nussinov. The function will
        (1) check if the letters of the input sequence are valid,
        (2) initialize the dynamic programming matrix self.matrix
        (3) calculate the optimal pairing
        (4) return the bracket representation of the optimal structure.
        :param sequence: input sequence
        :return: optimal structure.
        """
        self.alphabet.check_words(sequence)
        size = len(sequence)
        self.structure = []
        self.matrix = np.zeros((size, size))
        # fill the matrix
        for n in range(1, size):
            for j in range(n, size):
                i = j - n
                case1 = self.matrix[i + 1][j - 1] + 1 if self.is_basepair(sequence[i], sequence[j]) else 0
                case2 = self.matrix[i + 1][j]
                case3 = self.matrix[i][j - 1]
                if i + 3 <= j:
                    tmp = []
                    for k in range(i + 1, j):
                        tmp.append(self.matrix[i, k] + self.matrix[k + 1, j])
                    case4 = max(tmp)
                    self.matrix[i][j] = max(case1, case2, case3, case4)
                else:
                    self.matrix[i][j] = max(case1, case2, case3)
        # mirror
        for i in range(size):
            for j in range(0, i):
                self.matrix[i][j] = self.matrix[j][i]

        self.structure = self.traceback(sequence, 0, size - 1, [])
        structure, coords_list = self.structure_to_brackets_coords(sequence)
        return sequence, structure, coords_list


def process_program_arguments():
    if not args.input:
        LOGGER.critical(
                "Error, provide input file/files/directory/directories using -i / --input. -h/--help for all "
                "arguments.")
        exit(1)


def run_nussinov():
    sequences = parse_input(args.input, args.file_filter)
    nus = Nussinov(min_loop_length=args.min_loop_length)
    for sequence in sequences:
        sequence, brackets, coords_list = nus.run(sequence)
        LOGGER.info(f'Sequence: {sequence},\n'
                    f'Structure (Brackets): {brackets}\n'
                    f'Structure (coords): {coords_list}')


def main():
    process_program_arguments()
    run_nussinov()
    exit(1)


if __name__ == '__main__':
    # command line arguments.
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input', type=str, required=True, default='',
                        nargs='+',
                        help='A file, a directory or multiple directories. directories are processed recursively.')
    parser.add_argument('--file-filter', type=str, default='',
                        help='A regex to define a filter, which will be applied to the files parsed.')
    parser.add_argument('-l', '--min_loop_length', type=int, default=4,
                        help='min loop length')
    parser.add_argument('-v', '--verbose', action='store_true', default=False,
                        help='verbose output.')

    args = parser.parse_args()
    main()
