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

    def is_basepair(self, letter1, letter2):
        return {letter1, letter2} in self.basepairs

    def calc_optimal_pairing(self, i, j, sequence):
        """ returns the score of the optimal pairing between indices i and j"""
        # base case: no pairs allowed when i and j are less than 4 bases apart
        if i >= j - self.min_loop_length:
            return 0
        else:
            # i and j can be paired or not, if not paired then the optimal score is the optimal pairing at i,j-1
            not_paired = self.calc_optimal_pairing(i, j - 1, sequence)

            # check if j can be involved in a pairing with a position t
            paired = [1
                      + self.calc_optimal_pairing(i, k - 1, sequence)
                      + self.calc_optimal_pairing(k + 1, j - 1, sequence)
                      for k in range(i, j - 4) if self.is_basepair(sequence[k], sequence[j])
                      ]
            if not paired:
                paired = [0]
            paired = max(paired)
            return max(not_paired, paired)

    def traceback(self, i, j, sequence):
        """
        Traceback in the scoring matrix
        :param i: i pos in the matrix
        :param j: j pos in the matrix
        :param sequence: input sequence
        :return: void
        """
        # in this case we've gone through the whole sequence. Nothing to do.
        if j <= i:
            return
        # if j is unpaired, there will be no change in score when we take it out, so we just recurse to the next index
        elif self.matrix[i][j] == self.matrix[i][j - 1]:
            self.traceback(i, j - 1, sequence)
        else:
            # try pairing j with a matching index k to its left.
            for k in [b for b in range(i, j - self.min_loop_length) if self.is_basepair(sequence[b], sequence[j])]:
                # if the score at i,j is the result of adding 1 from pairing (j,k) and whatever score
                # comes from the substructure to its left (i, k-1) and to its right (k+1, j-1)
                if k - 1 < 0:
                    if self.matrix[i][j] == self.matrix[k + 1][j - 1] + 1:
                        self.structure.append((k, j))
                        self.traceback(k + 1, j - 1, sequence)
                elif self.matrix[i][j] == self.matrix[i][k - 1] + self.matrix[k + 1][j - 1] + 1:
                    # add the pair (j,k) to the list of pairs
                    self.structure.append((k, j))
                    # recursive exploration of substructures formed by this pairing
                    self.traceback(i, k - 1, sequence)
                    self.traceback(k + 1, j - 1, sequence)

    def structure_to_brackets(self, sequence):
        """
        Convert a structure to bracket representation
        :param sequence: input sequence.
        :return:  bracket representation
        """
        assert self.structure, "structure is empty."
        dot_list = ["." for _ in range(len(sequence))]
        for s in self.structure:
            dot_list[min(s)] = "("
            dot_list[max(s)] = ")"
        return "".join(dot_list)

    def structure_to_coords(self, sequence):
        """
        Convert a structure to bracket representation
        :param sequence: input sequence.
        :return:  bracket representation
        """
        assert self.structure, "structure is empty."
        coords = []
        for s in self.structure:
            coords.append((min(s), max(s)))
        return coords

    def init_matrix(self, size):
        """
        Initialize the matrix.
        :param size: lenght of the sequence.
        :return:
        """
        # size x size matrix.
        self.matrix = np.empty((size, size))
        self.matrix[:] = np.NAN
        for k in range(0, self.min_loop_length):
            for i in range(size - k):
                j = i + k
                self.matrix[i][j] = 0

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
        self.init_matrix(size=size)
        # fill the matrix
        for k in range(self.min_loop_length, size):
            for i in range(size - k):
                j = i + k
                self.matrix[i][j] = self.calc_optimal_pairing(i, j, sequence)

        # mirror
        for i in range(size):
            for j in range(0, i):
                self.matrix[i][j] = self.matrix[j][i]

        self.traceback(0, size - 1, sequence)
        return sequence, self.structure_to_brackets(sequence), self.structure_to_coords(sequence)


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
        sequence, brackets, coords = nus.run(sequence)
        LOGGER.info(f'Sequence: {sequence},\n'
                    f'Structure (Brackets): {brackets}\n'
                    f'Structure (coords): {coords}')


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
else:
    pass
