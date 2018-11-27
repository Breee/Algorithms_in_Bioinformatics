from enum import Enum

import numpy

from prakt.nw import NeedlemanWunschBase


class operation(Enum):
    MATCH = 1,
    INSERTION = 2,
    DELETION = 3,
    MISMATCH = 4,
    ROOT = 0


class TracebackCell(object):
    def __init__(self, type, predecessors, score):
        self.type = type
        self.predecessors = predecessors
        self.score = score

    def __str__(self):
        return "(%s, %s, %s)" % (self.type, self.predecessors, self.score)

    def __repr__(self):
        return "(%s, %s, %s)" % (self.type, self.predecessors, self.score)


class NeedlemanWunsch(NeedlemanWunschBase):
    """Document me!"""

    def __init__(self, match_scoring=1, indel_scoring=-1, mismatch_scoring=-1, substitution_matrix=None):
        self.match_scoring = match_scoring
        self.indel_scoring = indel_scoring
        self.mismatch_scoring = mismatch_scoring
        self.substitution_matrix = substitution_matrix
        self.traceback_matrix = numpy.zeros(shape=(1, 1), dtype=object)
        self.desired_traceback = TracebackCell(None, None, None)
        self.scoring_matrix = numpy.zeros(shape=(1, 1))
        self.tracebacks = dict()

    def calculate_distance_matrix(self, seq1, seq2):
        pass

    def init_scoring_matrix(self, seq1, seq2):
        # initialize scoring matrix with all zeros
        self.scoring_matrix = numpy.zeros(shape=(len(seq1), len(seq2)))
        self.traceback_matrix = numpy.zeros(shape=(len(seq1), len(seq2)), dtype=object)
        self.scoring_matrix[0][0] = 0
        self.traceback_matrix[0][0] = TracebackCell(type=operation.ROOT, predecessors=[], score=0)
        # iterate top row and initialize values
        for i in range(1, len(self.scoring_matrix[0])):
            score = self.mismatch_scoring * i
            self.scoring_matrix[0][i] = score
            self.traceback_matrix[0][i] = TracebackCell(type=operation.DELETION,
                                                        predecessors=[self.traceback_matrix[0][i - 1]], score=score)
            # iterate first column and initialize values
        for i in range(1, len(self.scoring_matrix.T[0])):
            score = self.mismatch_scoring * i
            self.scoring_matrix.T[0][i] = self.mismatch_scoring * i
            self.traceback_matrix.T[0][i] = TracebackCell(type=operation.INSERTION,
                                                          predecessors=[self.traceback_matrix.T[0][i - 1]], score=score)

    def score(self, letter1, letter2):
        if letter1 == letter2:
            return self.match_scoring
        else:
            return self.mismatch_scoring

    def calculate_similarity_matrix(self, seq1, seq2):
        """
        >>> nw = NeedlemanWunsch()
        >>> nw.calculate_similarity_matrix("€AATC","€AACT")
        pew

        :param seq1:
        :param seq2:
        :return:
        """
        # initialize scoring matrix.
        self.init_scoring_matrix(seq1=seq1, seq2=seq2)
        # next we want to fill the scoring matrix.
        for i in range(1, len(seq1)):
            for j in range(1, len(seq2)):
                # we calculate the score of the letters in the sequences (Match/Mismatch)
                score = self.score(seq1[i], seq2[j])
                # Top left cell
                match = self.scoring_matrix[i - 1][j - 1] + score
                # Top cell
                delete = self.scoring_matrix[i - 1][j] + self.indel_scoring
                # Left cell
                insert = self.scoring_matrix[i][j - 1] + self.indel_scoring

                # We calculate the maximum / minimum
                # If we want similarity, we maximize. If we want distance we minimize
                # TODO: Distance.
                maximum = max(match, delete, insert)
                self.scoring_matrix[i][j] = maximum
                cell = None
                if maximum == match:
                    cell = TracebackCell(type=operation.MATCH,
                                         predecessors=[self.traceback_matrix[i - 1][j - 1]],
                                         score=self.scoring_matrix[i][j])
                if maximum == delete:
                    tb_cell = TracebackCell(type=operation.DELETION,
                                            predecessors=[self.traceback_matrix[i - 1][j]],
                                            score=self.scoring_matrix[i][j])
                    if cell:
                        cell.predecessors.append(tb_cell)
                    else:
                        cell = tb_cell

                if maximum == insert:
                    tb_cell = TracebackCell(type=operation.INSERTION,
                                            predecessors=[self.traceback_matrix[i][j - 1]],
                                            score=self.scoring_matrix[i][j])
                    if cell:
                        cell.predecessors.append(tb_cell)
                    else:
                        cell = tb_cell

                self.traceback_matrix[i][j] = cell

        print(self.scoring_matrix)
        print(self.traceback_matrix[len(seq1) - 1][len(seq2) - 1])
        self.desired_traceback = self.traceback_matrix[len(seq1) - 1][len(seq2) - 1]
        self.split_traceback_set()
        print(self.tracebacks)
        self.generate_alignments(sequence1=seq1, sequence2=seq2)

    def explore(self, traceback_cell: TracebackCell, traceback_id):
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
        traceback_id = 0
        self.explore(self.desired_traceback, traceback_id)

    def generate_alignments(self, sequence1, sequence2):
        alignments = []
        for id, traceback in self.tracebacks.items():
            alignment = [[], []]
            i = 0
            j = 0
            for op in traceback:
                if op == operation.MATCH or op == operation.MISMATCH or op == operation.ROOT:
                    alignment[0].append(sequence1[i])
                    alignment[1].append(sequence2[j])
                    i += 1
                    j += 1
                elif op == operation.DELETION:
                    alignment[0].append(sequence1[i])
                    alignment[1].append("-")
                    i += 1
                elif op == operation.INSERTION:
                    alignment[0].append("-")
                    alignment[1].append(sequence2[j])
                    j += 1
            alignments.append(alignment)
        for alignment in alignments:
            print("ALIGNMENT:")
            print("%s\n%s" % ("".join(alignment[0]), "".join(alignment[1])))
            print("------------------")

    def run(self, seq1_fasta_fn, seq2_fasta_fn, subst_matrix_fn, cost_gap_open, complete_traceback):
        """Document me!"""

        # return some dummy results
        return ("idA",
                "FancySequenceA",
                "idB",
                "FancysequenceB",
                1000,
                [("Fancy_SequenceA_", "Fancys_equence_B")])


if __name__ == '__main__':
    # run Needleman-Wunsch with some parameters
    nw = NeedlemanWunsch()
    nw.run(
            "data/sequence1.fa",
            "data/sequence2.fa",
            "data/blosum62.txt",
            5,
            True)
