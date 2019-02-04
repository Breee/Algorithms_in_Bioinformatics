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
from pprint import pformat

import numpy
from Bio.SeqRecord import SeqRecord
from Bio.SubsMat import MatrixInfo

from logger.log import setup_custom_logger
from utility.utils import Alignment, Alphabet, Operation, Result, ScoringType, TracebackCell, parse_input

LOGGER = setup_custom_logger("nw", logfile="needleman_wunsch.log")

import argparse
import itertools


class ScoringSettings:
    def __init__(self, match_scoring=1, indel_scoring=-1, mismatch_scoring=-1, gap_penalty=6,
                 substitution_matrix=MatrixInfo.blosum62, similarity=True):
        self.similarity = similarity
        self.substitution_matrix = substitution_matrix
        self.gap_penalty = gap_penalty
        self.mismatch_scoring = mismatch_scoring
        self.indel_scoring = indel_scoring
        self.match_scoring = match_scoring


class NeedlemanWunsch(object):
    """
    Class which implements the needleman-wunsch algorithm.
    """

    def __init__(self, settings=ScoringSettings(), complete_traceback=False, verbose=False):
        LOGGER.info("Initialzing needleman-wunsch.")
        sigma = {"A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y",
                 "X"}
        self.alphabet = Alphabet(sigma)
        # scores
        self.match_scoring = settings.match_scoring
        self.indel_scoring = settings.indel_scoring
        self.mismatch_scoring = settings.mismatch_scoring
        # Subsitution matrices (PAM/BLOSSOM), default is blossom62.
        self.substitution_matrix = settings.substitution_matrix
        # The traceback matrix contains TracebackCell objects,
        self.traceback_matrix = numpy.zeros(shape=(1, 1), dtype=object)
        # Will store the TracebakCell in the bottom right of the Traceback Matrix.
        self.desired_traceback = TracebackCell(None, None)
        # All tracebacks
        self.tracebacks = list()
        self.complete_traceback = complete_traceback
        # The scoring matrix, which is used to calculate the optimal alignment scores.
        self.scoring_matrix = numpy.zeros(shape=(1, 1))
        if settings.similarity:
            self.scoring_type = ScoringType.SIMILARITY
            self.gap_penalty = abs(settings.gap_penalty) * (-1)
        else:
            self.scoring_type = ScoringType.DISTANCE
            self.gap_penalty = abs(settings.gap_penalty)
            raise NotImplementedError("DISTANCE is not implemented yet.")
        self.alignments = []
        if verbose:
            LOGGER.level = logging.DEBUG
        LOGGER.info(f'Needleman Wunsch initialized with: {["%s: %s" % item for item in vars(self).items()]}')

    def init_scoring_matrix(self, seq1, seq2):
        """
        This function
        :param seq1: DNA/RNA sequence
        :param seq2: DNA/RNA sequence
        :return:
        """
        LOGGER.info("Initializing Scoring Matrix.")
        # initialize scoring matrix with all zeros
        self.scoring_matrix = numpy.zeros(shape=(len(seq1) + 1, len(seq2) + 1))
        self.traceback_matrix = numpy.zeros(shape=(len(seq1) + 1, len(seq2) + 1), dtype=object)
        self.scoring_matrix[0][0] = 0
        self.traceback_matrix[0][0] = TracebackCell(predecessors=[], score=0)
        # iterate top row and initialize values
        for i in range(1, len(self.scoring_matrix[0])):
            score = self.gap_penalty * i
            self.scoring_matrix[0][i] = score
            self.traceback_matrix[0][i] = TracebackCell(
                    predecessors=[(Operation.DELETION, self.traceback_matrix[0][i - 1])], score=score)
            # iterate first column and initialize values
        for i in range(1, len(self.scoring_matrix.T[0])):
            score = self.gap_penalty * i
            self.scoring_matrix.T[0][i] = score
            self.traceback_matrix.T[0][i] = TracebackCell(
                    predecessors=[(Operation.INSERTION, self.traceback_matrix.T[0][i - 1])], score=score)
        assert self.scoring_matrix.shape == (len(seq1) + 1, len(seq2) + 1)
        assert self.traceback_matrix.shape == (len(seq1) + 1, len(seq2) + 1)

    def score(self, letter1, letter2):
        LOGGER.debug("Calculating score S(%s,%s)" % (letter1, letter2))
        # if the letter is a special gap char X from Multiple sequence alignment, return 0.
        if letter1 == "X" or letter2 == "X":
            return 0
        # if we use a substituion matrix, retrieve the score of a letter pair (l_1, l_2).
        elif self.substitution_matrix:
            pair = (letter1, letter2)
            if pair not in self.substitution_matrix:
                return self.substitution_matrix[(tuple(reversed(pair)))]
            elif pair in self.substitution_matrix:
                return self.substitution_matrix[pair]
            else:
                raise KeyError("%s and %s  do not appear as pair in substitution matrix: %s" % (
                    letter1, letter2, self.substitution_matrix))
        # Else use match / mismatch values as scores.
        elif letter1 == letter2:
            return self.match_scoring
        else:
            return self.mismatch_scoring

    def calculate_scoring_matrix(self, seq1, seq2):
        """
        Function which calculates the scoring matrix using needleman-wunsch.
        :param seq1: First sequence.
        :param seq2: Second sequence
        :return: void

        >>> nw = NeedlemanWunsch()
        >>> nw.calculate_scoring_matrix("AATC","AACT")
        >>> nw.scoring_matrix
        array([[  0.,  -6., -12., -18., -24.],
               [ -6.,   4.,  -2.,  -8., -14.],
               [-12.,  -2.,   8.,   2.,  -4.],
               [-18.,  -8.,   2.,   7.,   7.],
               [-24., -14.,  -4.,  11.,   6.]])
        >>> nw.traceback_matrix[-1][-1]
        ([(<Operation.MATCH: (1,)>, ([(<Operation.MATCH: (1,)>, ([(<Operation.MATCH: (1,)>, ([(<Operation.MATCH: (1,\
)>, ([], 0))], 4.0))], 8.0))], 7.0))], 6.0)
        """
        # initialize scoring matrix.
        self.init_scoring_matrix(seq1=seq1, seq2=seq2)
        LOGGER.debug("Calculating Score Matrix.")
        # next we want to fill the scoring matrix.
        for i in range(1, len(seq1) + 1):
            for j in range(1, len(seq2) + 1):
                letter1 = seq1[i - 1]
                letter2 = seq2[j - 1]
                # we calculate the score of the letters in the sequences (Match/Mismatch)
                score = self.score(letter1, letter2)
                # Top left cell
                match = self.scoring_matrix[i - 1][j - 1] + score
                # Top cell
                delete = self.scoring_matrix[i][j - 1] + self.gap_penalty
                # Left cell
                insert = self.scoring_matrix[i - 1][j] + self.gap_penalty
                # We calculate the maximum / minimum
                # If we want similarity, we maximize. If we want distance we minimize
                if self.scoring_type == ScoringType.DISTANCE:
                    optimum = min(match, delete, insert)
                else:
                    optimum = max(match, delete, insert)
                self.scoring_matrix[i][j] = optimum
                # We then update the traceback matrix.
                pre = []
                if optimum == match:
                    pre.append({"traceback_i": i - 1, "traceback_j": j - 1, "operation": Operation.MATCH})
                if optimum == delete:
                    pre.append({"traceback_i": i, "traceback_j": j - 1, "operation": Operation.DELETION})
                if optimum == insert:
                    pre.append({"traceback_i": i - 1, "traceback_j": j, "operation": Operation.INSERTION})
                cell = self.create_traceback_cell(predecessors=pre, score_i=i, score_j=j)
                self.traceback_matrix[i][j] = cell

    def create_traceback_cell(self, predecessors, score_i, score_j):
        """
        Helper function to create and update TracebackCells
        :param cell: TracebackCell.
        :param predecessors: predecessor states. e.g. [{"traceback_i": i - 1, "traceback_j": j, "operation":
        Operation.DELETION}]
        :param score_i: position i in scoring matrix.
        :param score_j: positon j in scoring matrix.
        :return:
        """

        tb_cell = TracebackCell(predecessors=[],
                                score=self.scoring_matrix[score_i][score_j])
        for pre in predecessors:
            tb_cell.predecessors.append((pre["operation"],
                                         self.traceback_matrix[pre["traceback_i"]][pre["traceback_j"]])
                                        )
        return tb_cell

    def explore(self, traceback_cell: TracebackCell, current_traceback=[]):
        """
        Function to explore a TracebackCell recursively to obtain all tracebacks to the root.
        :param traceback_cell: TracebackCell from which you want to start backtracking.
        :return: void, results are stored in member self.tracebacks (dictionary)
        """
        pre = traceback_cell.predecessors
        if not pre:
            return [current_traceback]

        tracebacks = []
        for op, cell in pre:
            # for each new predecessor, create a new traceback.
            tracebacks.extend(self.explore(cell, current_traceback=[op] + current_traceback))
        return tracebacks

    def split_traceback_set(self):
        """Function to kickoff exploration of a TracebackCell"""
        LOGGER.info("Kicking off exploration of Traceback.")
        self.tracebacks = self.explore(self.desired_traceback, [])

    def generate_alignments(self, seq1: SeqRecord, seq2: SeqRecord, all=False):
        """
        Function which creates all optimal alignments.
        :param sequence1:
        :param sequence2:
        :return:
        """
        sequence1 = seq1.seq
        sequence2 = seq2.seq
        LOGGER.info("Generating Alignments.")
        alignments = []
        for traceback in self.tracebacks:
            if not all and len(alignments) >= 1:
                self.alignments = alignments
                return
            LOGGER.debug("Length of tb: %d" % len(traceback))
            LOGGER.debug("tb: %s" % traceback)
            al_seq1 = ""
            al_seq2 = ""
            i = 0
            j = 0
            for op in traceback:
                if op == Operation.MATCH or op == Operation.MISMATCH:
                    al_seq1 += sequence1[i]
                    al_seq2 += sequence2[j]
                    i += 1
                    j += 1
                elif op == Operation.DELETION:
                    al_seq1 += "-"
                    al_seq2 += sequence2[j]
                    j += 1
                elif op == Operation.INSERTION:
                    al_seq1 += sequence1[i]
                    al_seq2 += "-"
                    i += 1
            assert len(al_seq1) == len(al_seq2)
            alignments.append(
                    Alignment(sequence1=SeqRecord(seq=al_seq1, id=seq1.id, name=seq1.name),
                              sequence2=SeqRecord(seq=al_seq2, id=seq2.id, name=seq2.name),
                              operations=[],
                              score=self.desired_traceback.score))

        self.alignments = alignments
        return

    def run(self, seq1, seq2, complete_traceback=False):
        """
        :param fasta_files:
        :param complete_traceback:
        :return:

        >>> nw = NeedlemanWunsch()
        >>> fasta = ["data/test1.fa", "data/test2.fa"]
        >>> from utility.utils import parse_fasta_files
        >>> sequences = parse_fasta_files(fasta)
        >>> res = nw.run(sequences[0],sequences[1], complete_traceback=False)
        >>> res.alignments
        [Alignment: (ID: S1
        Name: S1
        Description: <unknown description>
        Number of features: 0
        'AAA', ID: S2
        Name: S2
        Description: <unknown description>
        Number of features: 0
        '-AA'), Score: 2]
        """
        LOGGER.info("Running on sequences: (%s, %s)" % (seq1.seq, seq2.seq))
        self.alphabet.check_words({seq1.seq, seq2.seq})
        self.calculate_scoring_matrix(seq1.seq, seq2.seq)
        self.desired_traceback = self.traceback_matrix[-1][-1]
        self.split_traceback_set()
        self.generate_alignments(seq1=seq1, seq2=seq2, all=complete_traceback)

        res = None
        if self.alignments:
            if not complete_traceback:
                res = [self.alignments[0]]
            else:
                res = self.alignments

        # return some dummy results
        return Result(seq1_ID=seq1.id,
                      seq1=seq1,
                      seq2_ID=seq2.id,
                      seq2=seq2,
                      score=self.desired_traceback.score,
                      alignments=res)

    def pairwise_alignments(self, sequences):
        if len(sequences) < 2:
            LOGGER.warn("We received not enough sequences. Make sure you called the program correctly.")
            exit(1)
        results = []
        LOGGER.info("You provided a total of %d sequences, performing pairwise aligments." % len(sequences))
        combinations = list(itertools.combinations(sequences, 2))
        total = len(combinations)
        current = 1
        for x, y in sorted(combinations, key=lambda x: x[0].id):
            LOGGER.info(" Alignment %d / %d (SEQ1: %s, SEQ2: %s)" % (current, total, x.seq, y.seq))
            res = self.run(x, y, complete_traceback=self.complete_traceback)
            results.append(res)
            LOGGER.info("SEQUENCE PAIR:(ID:%s SEQ:%s) ~ (ID:%s SEQ:%s)" % (x.id, x.seq, y.id, y.seq))
            LOGGER.info("SCORE: %s" % res.score)
            LOGGER.info("PRINTING ALIGNMENT(S): ")
            for i, alignment in enumerate(res.alignments):
                LOGGER.info("%d. ALIGNMENT: (%s, %s)" % (i + 1, alignment.sequence1.seq, alignment.sequence2.seq))
            current += 1
        return results


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


def run_needleman():
    sequences = parse_input(args.input, args.file_filter)
    # init the needleman
    settings = ScoringSettings(substitution_matrix=args.substitution_matrix, gap_penalty=args.gap_penalty,
                               similarity=(not args.distance))
    nw = NeedlemanWunsch(settings, complete_traceback=args.all, verbose=args.verbose)
    results = nw.pairwise_alignments(sequences)
    LOGGER.info("SUMMARY:\n%s" % pformat(results))


def main():
    process_program_arguments()
    run_needleman()
    exit(1)


if __name__ == '__main__':
    # command line arguments.
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input', type=str, required=True, default='',
                        nargs='+',
                        help='A file, a directory or multiple directories. directories are processed recursively.')
    parser.add_argument('--file-filter', type=str, default='',
                        help='A regex to define a filter, which will be applied to the files parsed.')
    parser.add_argument('-a', '--all', action='store_true', default=False,
                        help='Return ALL optimal alignments.')
    parser.add_argument('-g', '--gap_penalty', type=float, default=6.0,
                        help='Gap penalty, float value. default is 6.0')
    parser.add_argument('-s', '--substitution_matrix', type=str, default="BLOSSUM",
                        help='Substitution Matrix (BLOSSUM | PAM | NONE) default is BLOSSUM')
    parser.add_argument('-d', '--distance', action='store_true', default=False,
                        help='Calculate DISTANCE instead of SIMILARITY')
    parser.add_argument('-v', '--verbose', action='store_true', default=False,
                        help='verbose output.')

    args = parser.parse_args()
    main()
