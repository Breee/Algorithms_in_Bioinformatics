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

import numpy

from logger.log import setup_custom_logger
from utility.utils import Alphabet, parse_input

LOGGER = setup_custom_logger("nov", logfile="nussinov.log")

import argparse


class Nussinov(object):
    """
    Class which implements the nussinov algorithm.
    """

    def __init__(self, min_loop_length=1, verbose=False):
        LOGGER.info("Initialzing nussinov")
        sigma = {"A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y",
                 "X", "U"}
        self.alphabet = Alphabet(sigma)
        self.matrix = None
        self.basepairs = {frozenset(["A", "U"]), frozenset(["G", "C"]), frozenset(["G", "U"])}
        self.min_loop_length = min_loop_length
        self.paired = dict()

    @staticmethod
    def kth_diag_indices(a, k):
        rows, cols = numpy.diag_indices_from(a)
        if k < 0:
            return rows[-k:], cols[:k]
        elif k > 0:
            return rows[:-k], cols[k:]
        else:
            return rows, cols

    def is_basepair(self, letter1, letter2):
        return {letter1, letter2} in self.basepairs

    def calculate_matrix(self, seq):
        """
        :param seq:
        :param loop_length:
        :return: void
        >>> nus = Nussinov()
        >>> nus.calculate_matrix("GCACGACG")
        >>> nus.matrix
        array([[0., 0., 1., 1., 1., 2., 2., 2., 3.],
               [0., 0., 0., 0., 0., 1., 1., 1., 2.],
               [0., 0., 0., 0., 0., 1., 1., 1., 2.],
               [0., 0., 0., 0., 0., 1., 1., 1., 2.],
               [0., 0., 0., 0., 0., 0., 0., 1., 1.],
               [0., 0., 0., 0., 0., 0., 0., 0., 1.],
               [0., 0., 0., 0., 0., 0., 0., 0., 1.],
               [0., 0., 0., 0., 0., 0., 0., 0., 0.]])
        """
        self.matrix = numpy.zeros(shape=(len(seq), len(seq) + 1), dtype=float)

        for i in range(2, len(seq) + 1):
            for k, j in zip(range(i, len(seq) + 1), range(0, len(seq) - 1)):
                self.calc_score(j, k, seq)

    def calc_score(self, i, j, seq):
        max_val = 0
        k = i
        while i <= k and k < j:
            if self.is_basepair(seq[k], seq[j - 1]):
                energy_of_pairing = self.matrix[i][k - 1] + self.matrix[k + 1][j - 1] + 1
                if max_val < energy_of_pairing:
                    max_val = energy_of_pairing
            k += 1
        self.matrix[i][j] = max(self.matrix[i][j - 1], max_val)

    def traceback(self, i, j, seq):
        if j <= i:
            return
        elif self.matrix[i][j] == self.matrix[i][j - 1]:
            self.traceback(i, j - 1, seq)
            return
        else:
            k = i
            while i <= k and k < j:
                if self.is_basepair(seq[k - 1], seq[j - 1]):
                    if self.matrix[i][j] == self.matrix[i][k - 1] + self.matrix[k][j - 1] + 1:
                        self.paired[k] = j
                        self.traceback(i, k - 1, seq)
                        self.traceback(k, j - 1, seq)
                        return
                k += 1

    def run(self, sequence):
        self.alphabet.check_words(sequence)
        self.calculate_matrix(sequence)
        self.traceback(0, len(sequence), sequence)
        print(self.paired)


def process_program_arguments():
    if not args.input:
        LOGGER.critical(
                "Error, provide input file/files/directory/directories using -i / --input. -h/--help for all "
                "arguments.")
        exit(1)


def run_nussinov():
    sequences = parse_input(args.input, args.file_filter)
    nus = Nussinov(min_loop_length=args.loop_length)
    for sequence in sequences:
        nus.run(sequence)


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
    parser.add_argument('-l', '--loop_length', type=int, default=1,
                        help='loop length')
    parser.add_argument('-v', '--verbose', action='store_true', default=False,
                        help='verbose output.')

    args = parser.parse_args()
    main()
else:
    pass
