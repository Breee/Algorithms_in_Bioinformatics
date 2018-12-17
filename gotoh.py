from pprint import pformat

import numpy
from Bio.SubsMat import MatrixInfo

from logger.log import setup_custom_logger
from utility.utils import Alignment, Operation, ScoringType, check_for_duplicates, parse_directory, parse_fasta_files, \
    split_directories_and_files

LOGGER = setup_custom_logger("nw", logfile="gotoh.log")

import argparse
import os
import itertools
import math


class TracebackCell(object):
    """
    A TracebackCell object which consists of
    predecessors: a list of TracebackCells, which are predecessors of this Cell.
    score:  score of this Cell.
    """

    def __init__(self, predecessors, score):
        self.predecessors = predecessors
        self.score = score

    def __str__(self):
        return "(%s, %s)" % (self.predecessors, self.score)

    def __repr__(self):
        return "(%s, %s)" % (self.predecessors, self.score)


class Result(object):
    """
    Represents a result tuple.
    """

    def __init__(self, seq1_ID, seq1, seq2_ID, seq2, alignments, score):
        self.seq1_ID = seq1_ID
        self.seq2_ID = seq2_ID
        self.seq1 = seq1
        self.seq2 = seq2
        self.alignments = alignments
        self.score = score

    def __repr__(self):
        return "(SEQ1: %s, %s, SEQ2: %s, %s,\n ALIGNMENTS:\n%s,\n SCORE: %s)" % (
            self.seq1_ID, self.seq1, self.seq2_ID, self.seq2, pformat(self.alignments), self.score)


class Gotoh(object):
    """
    Class which implements the needleman-wunsch algorithm.
    """

    def __init__(self, match_scoring=1, indel_scoring=-1, mismatch_scoring=-1, gap_penalty=-6, gap_extend=-3,
                 substitution_matrix=MatrixInfo.blosum62,
                 similarity=True):
        LOGGER.info("Initialzing needleman-wunsch.")
        # scores
        self.match_scoring = match_scoring
        self.indel_scoring = indel_scoring
        self.mismatch_scoring = mismatch_scoring
        # Subsitution matrices (PAM/BLOSSOM), default is blossom62.
        self.substitution_matrix = substitution_matrix
        # The traceback matrix contains TracebackCell objects,
        self.traceback_matrix = numpy.zeros(shape=(1, 1), dtype=object)
        # Will store the TracebakCell in the bottom right of the Traceback Matrix.
        self.desired_traceback = TracebackCell(None, None)
        # All tracebacks
        self.tracebacks = list()
        # The scoring matrix, which is used to calculate the optimal alignment scores.
        self.scoring_matrix_D = numpy.zeros(shape=(1, 1))
        self.scoring_matrix_Q = numpy.zeros(shape=(1, 1))
        self.scoring_matrix_P = numpy.zeros(shape=(1, 1))

        if similarity:
            self.scoring_type = ScoringType.SIMILARITY
            self.gap_penalty = gap_penalty
            self.gap_extend = gap_extend
        else:
            self.scoring_type = ScoringType.DISTANCE
            self.gap_penalty = abs(gap_penalty)
            self.gap_extend = abs(gap_extend)
            raise NotImplementedError("DISTANCE is not implemented yet.")
        self.alignments = []
        LOGGER.info("Scoring-Type: %s" % self.scoring_type)
        LOGGER.debug("Substitution Matrix:\n %s" % pformat(self.substitution_matrix))
        LOGGER.info("Gap Penalty: %d" % self.gap_penalty)

    def init_scoring_matrices(self, seq1, seq2):
        """
        This function
        :param seq1: DNA/RNA sequence
        :param seq2: DNA/RNA sequence
        :return:

        >>> got = Gotoh()
        >>> got.init_scoring_matrices("AAA","AAC")
        >>> got.scoring_matrix_D
        array([[  0.,  -9., -12., -15.],
               [ -9.,   0.,   0.,   0.],
               [-12.,   0.,   0.,   0.],
               [-15.,   0.,   0.,   0.]])
        >>> got.scoring_matrix_Q
        array([[  0., -inf, -inf, -inf],
               [ -9.,   0.,   0.,   0.],
               [-12.,   0.,   0.,   0.],
               [-15.,   0.,   0.,   0.]])
        >>> got.scoring_matrix_P
        array([[  0.,  -9., -12., -15.],
               [-inf,   0.,   0.,   0.],
               [-inf,   0.,   0.,   0.],
               [-inf,   0.,   0.,   0.]])
        """
        LOGGER.info("Initializing Scoring Matrix.")
        # initialize scoring matrix with all zeros
        self.scoring_matrix_D = numpy.zeros(shape=(len(seq1) + 1, len(seq2) + 1))
        self.scoring_matrix_Q = numpy.zeros(shape=(len(seq1) + 1, len(seq2) + 1))
        self.scoring_matrix_P = numpy.zeros(shape=(len(seq1) + 1, len(seq2) + 1))
        self.traceback_matrix = numpy.zeros(shape=(len(seq1) + 1, len(seq2) + 1), dtype=object)
        # next, set initial values according to gotohs algorithm.
        self.scoring_matrix_D[0][0] = 0
        # Q[0][0] is never used, set to it to -inf
        self.scoring_matrix_Q[0][0] = -math.inf
        # P[0][0] is never used, set to it to -inf
        self.scoring_matrix_P[0][0] = -math.inf
        # TB[0][0] is the root, thus no predecessors.
        self.traceback_matrix[0][0] = TracebackCell(predecessors=[], score=0)
        # next, iterate top row and initialize values
        for j in range(1, len(self.scoring_matrix_D[0])):
            score = self.gap_cost(j)
            self.scoring_matrix_D[0][j] = self.gap_cost(j)
            # P[0][j] and Q[0][j] is never used, set to it to -inf
            self.scoring_matrix_Q[0][j] = -math.inf
            self.scoring_matrix_P[0][j] = -math.inf
            self.traceback_matrix[0][j] = TracebackCell(
                    predecessors=[(Operation.DELETION, self.traceback_matrix[0][j - 1])], score=score)
            # iterate first column and initialize values
        for i in range(1, len(self.scoring_matrix_D.T[0])):
            score = self.gap_cost(i)
            self.scoring_matrix_D[i][0] = self.gap_cost(i)
            # P[i][0] and Q[i][0] is never used, set to it to -inf
            self.scoring_matrix_Q[i][0] = -math.inf
            self.scoring_matrix_P[i][0] = -math.inf
            self.traceback_matrix[i][0] = TracebackCell(
                    predecessors=[(Operation.INSERTION, self.traceback_matrix[i - 1][0])], score=score)
        LOGGER.info("Done.")

    def gap_cost(self, length):
        """
        affine gap cost function.
        :param length:
        :return: cost for gap of size :length:
        """
        return self.gap_penalty + (self.gap_extend * length)

    def score(self, letter1, letter2):
        """
        Scoring function, takes the score s(letter1,letter2) from self.substitution_matrix, if it is not None.
        If the self.substitution_matrix is None, then self.match_scoring + self.mismatch_scoring are used.

        :param letter1:
        :param letter2:
        :return: numeric score s(letter1,letter2)
        """
        LOGGER.debug("Calculating score S(%s,%s)" % (letter1, letter2))
        if self.substitution_matrix:
            pair = (letter1, letter2)
            if pair not in self.substitution_matrix:
                return self.substitution_matrix[(tuple(reversed(pair)))]
            elif pair in self.substitution_matrix:
                return self.substitution_matrix[pair]
            else:
                raise KeyError("%s and %s  do not appear as pair in substitution matrix: %s" % (
                    letter1, letter2, self.substitution_matrix))
        elif letter1 == letter2:
            return self.match_scoring
        else:
            return self.mismatch_scoring

    def calculate_scoring_matrix(self, seq1, seq2):
        """
        Function which calculates the scoring matrix using gotoh.
        Using the following recursion:

        1. P[i][j] = max(D[i-1][j] + g(1), P[i-1][j] + \beta)
        2. Q[i][j] = max(D[i][j-1] + g(1), Q[i][j-1] + \beta)
        3. D[i][j] = max(D[i-1][j-1] + score(x,y), Q[i][j], P[i][j])

        where:
        - score(x,y) is a scoring function,
        - \beta is the cost of a gap_extension
        - g(k) is the gap cost for a gap of size k

        :param seq1: First sequence.
        :param seq2: Second sequence
        :return: void. the function fills self.scoring_matrix_D, self.scoring_matrix_P, self.scoring_matrix_Q,
        self.traceback_matrix

        >>> got = Gotoh(match_scoring=1, mismatch_scoring=-1, substitution_matrix=None)
        >>> got.calculate_scoring_matrix("AATC","AACT")
        >>> got.scoring_matrix_D
        array([[  0.,  -9., -12., -15., -18.],
               [ -9.,   1.,  -8., -11., -14.],
               [-12.,  -8.,   2.,  -7., -10.],
               [-15., -11.,  -7.,   1.,  -6.],
               [-18., -14., -10.,  -6.,   0.]])
        >>> got.scoring_matrix_P
        array([[-inf, -inf, -inf, -inf, -inf],
               [-inf, -18., -21., -24., -27.],
               [-inf,  -8., -17., -20., -23.],
               [-inf, -11.,  -7., -16., -19.],
               [-inf, -14., -10.,  -8., -15.]])
        >>> got.scoring_matrix_Q
        array([[-inf, -inf, -inf, -inf, -inf],
               [-inf, -18.,  -8., -11., -14.],
               [-inf, -21., -17.,  -7., -10.],
               [-inf, -24., -20., -16.,  -8.],
               [-inf, -27., -23., -19., -15.]])
        """
        # initialize scoring matrix.
        self.init_scoring_matrices(seq1=seq1, seq2=seq2)
        LOGGER.debug("Calculating Score Matrix.")
        # next we want to fill the scoring matrix.
        for i in range(1, len(seq1) + 1):
            for j in range(1, len(seq2) + 1):
                letter1 = seq1[i - 1]
                letter2 = seq2[j - 1]
                # we calculate the score of the letters in the sequences (Match/Mismatch)
                score = self.score(letter1, letter2)

                # calculate Pi,j
                # calculate Qi,j
                # calculate Di,j
                self.scoring_matrix_P[i][j] = max(self.scoring_matrix_D[i - 1][j] + self.gap_cost(1),
                                                  self.scoring_matrix_P[i - 1][j] + self.gap_extend)
                self.scoring_matrix_Q[i][j] = max(self.scoring_matrix_D[i][j - 1] + self.gap_cost(1),
                                                  self.scoring_matrix_Q[i][j - 1] + self.gap_extend)

                match_D = self.scoring_matrix_D[i - 1][j - 1] + score
                insertion_Q = self.scoring_matrix_Q[i][j]
                deletion_P = self.scoring_matrix_P[i][j]
                self.scoring_matrix_D[i][j] = max(match_D, insertion_Q, deletion_P)

                # We then update the traceback matrix.
                pre = []
                if self.scoring_matrix_D[i][j] == match_D:
                    pre.append({"traceback_i": i - 1, "traceback_j": j - 1, "operation": Operation.MATCH})
                if self.scoring_matrix_D[i][j] == deletion_P:
                    pre.append({"traceback_i": i - 1, "traceback_j": j, "operation": Operation.DELETION})
                if self.scoring_matrix_D[i][j] == insertion_Q:
                    pre.append({"traceback_i": i, "traceback_j": j - 1, "operation": Operation.INSERTION})
                cell = self.create_traceback_cell(predecessors=pre, score_i=i, score_j=j)
                self.traceback_matrix[i][j] = cell

    def create_traceback_cell(self, predecessors, score_i, score_j):
        """
        Helper function to create and TracebackCells and fill them with predecessors.

        :param cell: TracebackCell.
        :param predecessors: predecessor states. e.g. [{"traceback_i": i - 1, "traceback_j": j, "operation":
        Operation.DELETION}]
        :param score_i: position i in scoring matrix.
        :param score_j: positon j in scoring matrix.
        :return: TracebackCell object.
        """

        tb_cell = TracebackCell(predecessors=[],
                                score=self.scoring_matrix_D[score_i][score_j])
        for pre in predecessors:
            tb_cell.predecessors.append((pre["operation"],
                                         self.traceback_matrix[pre["traceback_i"]][pre["traceback_j"]])
                                        )
        return tb_cell

    def explore(self, traceback_cell: TracebackCell, current_traceback=list()):
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
        LOGGER.info("Done.")

    def generate_alignments(self, sequence1, sequence2, all=False):
        """
        Function which creates all optimal alignments.
        :param sequence1:
        :param sequence2:
        :return:
        """
        LOGGER.info("Generating Alignments.")
        alignments = []
        for traceback in self.tracebacks:
            if not all and len(alignments) >= 1:
                self.alignments = alignments
                return
            LOGGER.debug("Length of tb: %d" % len(traceback))
            LOGGER.debug("tb: %s" % traceback)
            seq1 = ""
            seq2 = ""
            i = 0
            j = 0
            for op in traceback:
                if op == Operation.MATCH or op == Operation.MISMATCH:
                    seq1 += sequence1[i]
                    seq2 += sequence2[j]
                    i += 1
                    j += 1
                elif op == Operation.INSERTION:
                    seq1 += "-"
                    seq2 += sequence2[j]
                    j += 1
                elif op == Operation.DELETION:
                    seq1 += sequence1[i]
                    seq2 += "-"
                    i += 1
            alignments.append(Alignment(sequence1=seq1, sequence2=seq2,
                                        operations=[],
                                        score=self.desired_traceback.score))
        self.alignments = alignments
        return

    def run(self, seq1, seq2, complete_traceback=False):
        """
        :param fasta_files:
        :param complete_traceback:
        :return:

        >>> got = Gotoh(match_scoring=1, mismatch_scoring=-1,substitution_matrix=None)
        >>> fasta = ["data/test1.fa", "data/test2.fa"]
        >>> sequences = parse_fasta_files(fasta)
        >>> res = got.run(sequences[0],sequences[1], complete_traceback=True)
        >>> res
        (SEQ1: test1, ACA, SEQ2: test2, AA,
         ALIGNMENTS:
        [Alignment: (ACA, A-A), Score: -7],
         SCORE: -7.0)
        """
        self.calculate_scoring_matrix(seq1.seq, seq2.seq)
        self.desired_traceback = self.traceback_matrix[-1][-1]
        self.split_traceback_set()
        self.generate_alignments(sequence1=seq1.seq, sequence2=seq2.seq, all=complete_traceback)

        res = None
        if self.alignments:
            if not complete_traceback:
                res = [self.alignments[0]]
            else:
                res = self.alignments

        # return some dummy results
        return Result(seq1_ID=seq1.id,
                      seq1=seq1.seq,
                      seq2_ID=seq2.id,
                      seq2=seq2.seq,
                      score=self.desired_traceback.score,
                      alignments=res)

    def pairwise_alignments(self, sequences):
        results = []
        LOGGER.info("You provided a total of %d sequences, performing pairwise aligments." % len(sequences))
        combinations = list(itertools.combinations(sequences, 2))
        total = len(combinations)
        current = 1
        for x, y in combinations:
            LOGGER.info(" Alignment %d / %d (SEQ1: %s, SEQ2: %s)" % (current, total, x.seq, y.seq))
            res = self.run(x, y, complete_traceback=args.all)
            results.append(res)
            LOGGER.info("SEQUENCE PAIR:\nSEQ1 ID:%s SEQ:%s\n"
                        "SEQ2 ID:%s SEQ:%s" % (x.id, x.seq, y.id, y.seq))
            LOGGER.info("SCORE: %s" % res.score)
            LOGGER.info("PRINTING ALIGNMENT(S): ")
            for alignment in res.alignments:
                LOGGER.info("ALIGNMENT\n%s\n%s" % (alignment.sequence1, alignment.sequence2))
                print()
            current += 1
        return results


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
    if args.substitution_matrix != "PAM" and args.substitution_matrix != "BLOSSUM":
        LOGGER.critical(
                "UNKNOWN parameter for substitution matrix. Choose BLOSSUM or PAM, falling back to BLOSSUM")
        args.substitution_matrix = MatrixInfo.blosum62
    elif args.substitution_matrix == "BLOSSUM":
        args.substitution_matrix = MatrixInfo.blosum62
    elif args.substitution_matrix == "PAM":
        args.substitution_matrix = MatrixInfo.pam250


def run_gotoh():
    sequences = parse_input()
    if len(sequences) < 2:
        LOGGER.warn("We received not enough sequences. Make sure you called the program correctly.")
        exit(1)
    elif len(sequences) >= 2:
        # init the gotoh
        nw = Gotoh(substitution_matrix=args.substitution_matrix, gap_penalty=args.gap_penalty,
                   similarity=(not args.distance))
        results = nw.pairwise_alignments(sequences)
        LOGGER.info("SUMMARY:\n%s" % pformat(results))


def main():
    process_program_arguments()
    run_gotoh()
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
                        help='Substitution Matrix (BLOSSUM | PAM) default is BLOSSUM')
    parser.add_argument('-d', '--distance', action='store_true', default=False,
                        help='Calculate DISTANCE instead of SIMILARITY')

    args = parser.parse_args()
    main()
else:
    pass
