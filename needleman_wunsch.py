from pprint import pformat

import numpy
from Bio.SubsMat import MatrixInfo

from logger.log import setup_custom_logger
from utility.utils import Alignment, Alphabet, Operation, ScoringType, check_for_duplicates, parse_directory, \
    parse_fasta_files, split_directories_and_files

LOGGER = setup_custom_logger("nw", logfile="needleman_wunsch.log")

import argparse
import os
import itertools


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


class NeedlemanWunsch(object):
    """
    Class which implements the needleman-wunsch algorithm.
    """

    def __init__(self, match_scoring=1, indel_scoring=-1, mismatch_scoring=-1, gap_penalty=6,
                 substitution_matrix=MatrixInfo.blosum62,
                 similarity=True):
        LOGGER.info("Initialzing needleman-wunsch.")
        sigma = {"A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y"}
        self.alphabet = Alphabet(sigma)
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
        self.scoring_matrix = numpy.zeros(shape=(1, 1))
        if similarity:
            self.scoring_type = ScoringType.SIMILARITY
            self.gap_penalty = abs(gap_penalty) * (-1)
        else:
            self.scoring_type = ScoringType.DISTANCE
            self.gap_penalty = abs(gap_penalty)
            raise NotImplementedError("DISTANCE is not implemented yet.")
        self.alignments = []
        LOGGER.info("Scoring-Type: %s" % self.scoring_type)
        LOGGER.debug("Substitution Matrix:\n %s" % pformat(self.substitution_matrix))
        LOGGER.info("Gap Penalty: %d" % self.gap_penalty)

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
        LOGGER.info("Done.")
        assert self.scoring_matrix.shape == (len(seq1) + 1, len(seq2) + 1)
        assert self.traceback_matrix.shape == (len(seq1) + 1, len(seq2) + 1)

    def score(self, letter1, letter2):
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
                elif op == Operation.DELETION:
                    seq1 += "-"
                    seq2 += sequence2[j]
                    j += 1
                elif op == Operation.INSERTION:
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

        >>> nw = NeedlemanWunsch()
        >>> fasta = ["data/test1.fa", "data/test2.fa"]
        >>> sequences = parse_fasta_files(fasta)
        >>> res = nw.run(sequences[0],sequences[1], complete_traceback=False)
        >>> res
        ('test1', Seq('AAAAA', SingleLetterAlphabet()), 'test2', Seq('AAAA', SingleLetterAlphabet()), 10.0, \
[Alignment: (AAAAA, -AAAA), Score: 10])
        """
        self.alphabet.check_words({seq1, seq2})
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


def run_needleman():
    sequences = parse_input()
    if len(sequences) < 2:
        LOGGER.warn("We received not enough sequences. Make sure you called the program correctly.")
        exit(1)
    elif len(sequences) >= 2:
        # init the needleman
        nw = NeedlemanWunsch(substitution_matrix=args.substitution_matrix, gap_penalty=args.gap_penalty,
                             similarity=(not args.distance))
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
                        help='Substitution Matrix (BLOSSUM | PAM) default is BLOSSUM')
    parser.add_argument('-d', '--distance', action='store_true', default=False,
                        help='Calculate DISTANCE instead of SIMILARITY')

    args = parser.parse_args()
    main()
else:
    pass
