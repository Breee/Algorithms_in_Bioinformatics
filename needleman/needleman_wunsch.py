import logging
from pprint import pformat

import numpy
from Bio.SubsMat import MatrixInfo

from logger import log
from utility.utils import Alignment, Operation, ScoringType, parse_fasta_files

LOGGER = log.setup_custom_logger("needleman_wunsch", logfile="logfile.log")


class TracebackCell(object):
    """
    A TracebackCell object which consists of
    type: Type of Operation (MATCH/INSERTION/DELETION/MISMATCH/ROOT) of this Cell.
    predecessors: a list of TracebackCells, which are predecessors of this Cell.
    score:  score of this Cell.
    """

    def __init__(self, type, predecessors, score):
        self.type = type
        self.predecessors = predecessors
        self.score = score

    def __str__(self):
        return "(%s, %s, %s)" % (self.type, self.predecessors, self.score)

    def __repr__(self):
        return "(%s, %s, %s)" % (self.type, self.predecessors, self.score)


class NeedlemanWunsch(object):
    """
    Class which implements the needleman-wunsch algorithm.
    """

    def __init__(self, match_scoring=1, indel_scoring=-1, mismatch_scoring=-1, gap_penalty=6,
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
        self.desired_traceback = TracebackCell(None, None, None)
        # All tracebacks
        self.tracebacks = dict()
        # The scoring matrix, which is used to calculate the optimal alignment scores.
        self.scoring_matrix = numpy.zeros(shape=(1, 1))
        if similarity:
            self.scoring_type = ScoringType.SIMILARITY
            self.gap_penalty = abs(gap_penalty) * (-1)
        else:
            self.scoring_type = ScoringType.DISTANCE
            self.gap_penalty = abs(gap_penalty)
        self.alignments = []
        LOGGER.info("Scoring-Type: %s" % self.scoring_type)
        LOGGER.info("Substitution Matrix:\n %s" % pformat(self.substitution_matrix))
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
        self.scoring_matrix = numpy.zeros(shape=(len(seq1), len(seq2)))
        self.traceback_matrix = numpy.zeros(shape=(len(seq1), len(seq2)), dtype=object)
        self.scoring_matrix[0][0] = 0
        self.traceback_matrix[0][0] = TracebackCell(type=Operation.ROOT, predecessors=[], score=0)
        # iterate top row and initialize values
        for i in range(1, len(self.scoring_matrix[0])):
            score = self.gap_penalty * i
            self.scoring_matrix[0][i] = score
            self.traceback_matrix[0][i] = TracebackCell(type=Operation.DELETION,
                                                        predecessors=[self.traceback_matrix[0][i - 1]], score=score)
            # iterate first column and initialize values
        for i in range(1, len(self.scoring_matrix.T[0])):
            score = self.gap_penalty * i
            self.scoring_matrix.T[0][i] = score
            self.traceback_matrix.T[0][i] = TracebackCell(type=Operation.INSERTION,
                                                          predecessors=[self.traceback_matrix.T[0][i - 1]], score=score)
        LOGGER.info("Done.")
        assert self.scoring_matrix.shape == (len(seq1), len(seq2))
        assert self.traceback_matrix.shape == (len(seq1), len(seq2))

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
        >>> nw = NeedlemanWunsch()
        >>> nw.calculate_scoring_matrix("AATC","AACT")
        >>> nw.scoring_matrix
        array([[  0.,  -6., -12., -18.],
               [ -6.,   4.,   3.,   2.],
               [-12.,   3.,   3.,   8.],
               [-18.,   2.,  12.,  11.]])


        Function which calculates the scoring matrix using needleman-wunsch.
        :param seq1: First sequence.
        :param seq2: Second sequence
        :return: void
        """
        # initialize scoring matrix.
        self.init_scoring_matrix(seq1=seq1, seq2=seq2)
        LOGGER.info("Calculating Score Matrix.")
        # next we want to fill the scoring matrix.
        for i in range(1, len(seq1)):
            for j in range(1, len(seq2)):
                letter1 = seq1[i]
                letter2 = seq2[j]
                # we calculate the score of the letters in the sequences (Match/Mismatch)
                score = self.score(letter1, letter2)
                # Top left cell
                match = self.scoring_matrix[i - 1][j - 1] + score
                # Top cell
                delete = self.scoring_matrix[i - 1][j] + self.gap_penalty
                # Left cell
                insert = self.scoring_matrix[i][j - 1] + self.gap_penalty
                # We calculate the maximum / minimum
                # If we want similarity, we maximize. If we want distance we minimize
                if self.scoring_type == ScoringType.DISTANCE:
                    optimum = min(match, delete, insert)
                else:
                    optimum = max(match, delete, insert)
                self.scoring_matrix[i][j] = optimum
                # We then update the traceback matrix.
                cell = None
                pre = []
                if optimum == match:
                    pre.append((i - 1, j - 1, Operation.MATCH))
                    cell = self.create_traceback_cell(cell=cell, operation=Operation.MATCH, traceback_i=i - 1,
                                                      traceback_j=j - 1, score_i=i, score_j=j)
                if optimum == delete:
                    pre.append((i - 1, j, Operation.DELETION))
                    cell = self.create_traceback_cell(cell=cell, operation=Operation.DELETION, traceback_i=i - 1,
                                                      traceback_j=j, score_i=i, score_j=j)
                if optimum == insert:
                    pre.append((i, j - 1, Operation.INSERTION))
                    cell = self.create_traceback_cell(cell=cell, operation=Operation.INSERTION, traceback_i=i,
                                                      traceback_j=j - 1, score_i=i, score_j=j)
                print(pre)
                self.traceback_matrix[i][j] = cell
        print(self.scoring_matrix)

    def create_traceback_cell(self, cell, operation, traceback_i, traceback_j, score_i, score_j):
        """
        Helper function to create and update TracebackCells
        :param cell: TracebackCell.
        :param operation: Type of operation (Operation enum).
        :param traceback_i: position i in traceback_matrix of predecessor.
        :param traceback_j: position j in traceback_matrix of predecessor.
        :param score_i: position i in scoring matrix.
        :param score_j: positon j in scoring matrix.
        :return:
        """
        tb_cell = TracebackCell(type=operation,
                                predecessors=[self.traceback_matrix[traceback_i][traceback_j]],
                                score=self.scoring_matrix[score_i][score_j])
        if cell:
            cell.predecessors.append(tb_cell)
        else:
            cell = tb_cell
        return cell

    def explore(self, traceback_cell: TracebackCell, traceback_id=0):
        """
        Function to explore a TracebackCell recursively to obtain all tracebacks to the root.
        :param traceback_cell: TracebackCell from which you want to start backtracking.
        :param traceback_id: ID of the current traceback you explore.
        :return: void, results are stored in member self.tracebacks (dictionary)
        """
        if traceback_id in self.tracebacks:
            self.tracebacks[traceback_id].insert(0, traceback_cell.type)
        else:
            self.tracebacks[traceback_id] = [traceback_cell.type]
        pre = traceback_cell.predecessors
        if len(pre) > 1:
            for cell in pre:
                self.explore(cell, traceback_id)
                traceback_id += 1
        elif pre:
            self.explore(pre[0], traceback_id)

    def split_traceback_set(self):
        """Function to kickoff exploration of a TracebackCell"""
        LOGGER.info("Kicking off exploration of Traceback.")
        traceback_id = 0
        self.explore(self.desired_traceback, traceback_id)
        LOGGER.info("Done.")

    def generate_alignments(self, sequence1, sequence2):
        """
        Function which creates all optimal alignments.
        :param sequence1:
        :param sequence2:
        :return:
        """
        alignments = []
        for id, traceback in self.tracebacks.items():
            LOGGER.info("Length of tb: %d" % len(traceback))
            alignment = [[], [], []]
            i = 0
            j = 0
            for op in traceback:
                if op == Operation.MATCH or op == Operation.MISMATCH or op == Operation.ROOT:
                    alignment[0].append(sequence1[i])
                    alignment[1].append(sequence2[j])
                    i += 1
                    j += 1
                elif op == Operation.DELETION:
                    alignment[0].append(sequence1[i])
                    alignment[1].append("-")
                    i += 1
                elif op == Operation.INSERTION:
                    alignment[0].append("-")
                    alignment[1].append(sequence2[j])
                    j += 1
                alignment[2].append(op)
                if LOGGER.level == logging.DEBUG:
                    LOGGER.debug("######################")
                    LOGGER.debug("op: %s" % op)
                    LOGGER.debug("i:%d" % i)
                    LOGGER.debug("j:%d" % j)
                    LOGGER.debug("".join(alignment[0]), "len: %d" % len(alignment[0]))
                    LOGGER.debug("".join(alignment[1]), "len: %d" % len(alignment[1]))
                    LOGGER.debug("######################")
            alignments.append(
                    Alignment(sequence1="".join(alignment[0][1:]), sequence2="".join(alignment[1][1:]),
                              operations=alignment[2][1:],
                              score=self.desired_traceback.score))
        self.alignments = alignments

    def run(self, seq1, seq2, complete_traceback=False):
        """
        >>> nw = NeedlemanWunsch()
        >>> fasta = ["../data/test1.fn", "../data/test2.fn"]
        >>> sequences = parse_fasta_files(fasta)
        >>> res = nw.run(sequences[0],sequences[1])
        >>> res
        ('test1', Seq('€AAAC', SingleLetterAlphabet()), 'test2', Seq('€AAAG', SingleLetterAlphabet()), 10.0, (AAA-C, AAAG-, 10))


        :param fasta_files:
        :param complete_traceback:
        :return:
        """
        # add character wih represents the empty word
        seq1.seq = "€" + seq1.seq
        seq2.seq = "€" + seq2.seq
        self.calculate_scoring_matrix(seq1.seq, seq2.seq)
        self.desired_traceback = self.traceback_matrix[len(seq1.seq) - 1][len(seq2.seq) - 1]
        self.split_traceback_set()
        self.generate_alignments(sequence1=seq1.seq, sequence2=seq2.seq)

        if not complete_traceback:
            res = self.alignments[0]
        else:
            res = self.alignments

        # return some dummy results
        return (seq1.id,
                seq1.seq,
                seq2.id,
                seq2.seq,
                self.desired_traceback.score,
                res)


if __name__ == '__main__':
    blossum62 = MatrixInfo.blosum62
    # run Needleman-Wunsch with some parameters
    nw = NeedlemanWunsch()
    fasta_files = ["../data/test3.fn"]
    sequences = parse_fasta_files(fasta_files)
    seq1 = sequences[0]
    seq2 = sequences[1]
    res = nw.run(seq1, seq2, complete_traceback=True)
    LOGGER.debug("scoring:", nw.scoring_matrix)
    LOGGER.debug("Traceback-matrix:\n %s" % pformat(nw.traceback_matrix))
    LOGGER.debug("Traceback:\n %s" % pformat(nw.desired_traceback))
    LOGGER.debug("Done.")
    print(res)
