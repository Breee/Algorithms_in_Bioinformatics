import numpy
from Bio.SubsMat import MatrixInfo

from utility.utils import Alignment, Operation, ScoringType, parse_fasta_files


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
        # scores
        self.match_scoring = match_scoring
        self.indel_scoring = indel_scoring
        self.mismatch_scoring = mismatch_scoring
        self.gap_penalty = gap_penalty
        # Subsitution matrices (PAM/BLOSSOM).
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
        else:
            self.scoring_type = ScoringType.DISTANCE
        self.alignments = []

    def init_scoring_matrix(self, seq1, seq2):
        """
        This function
        :param seq1: DNA/RNA sequence
        :param seq2: DNA/RNA sequence
        :return:
        """
        # initialize scoring matrix with all zeros
        self.scoring_matrix = numpy.zeros(shape=(len(seq1), len(seq2)))
        self.traceback_matrix = numpy.zeros(shape=(len(seq1), len(seq2)), dtype=object)
        self.scoring_matrix[0][0] = 0
        self.traceback_matrix[0][0] = TracebackCell(type=Operation.ROOT, predecessors=[], score=0)
        # iterate top row and initialize values
        for i in range(1, len(self.scoring_matrix[0])):
            score = self.mismatch_scoring * i
            self.scoring_matrix[0][i] = score
            self.traceback_matrix[0][i] = TracebackCell(type=Operation.DELETION,
                                                        predecessors=[self.traceback_matrix[0][i - 1]], score=score)
            # iterate first column and initialize values
        for i in range(1, len(self.scoring_matrix.T[0])):
            score = self.mismatch_scoring * i
            self.scoring_matrix.T[0][i] = self.mismatch_scoring * i
            self.traceback_matrix.T[0][i] = TracebackCell(type=Operation.INSERTION,
                                                          predecessors=[self.traceback_matrix.T[0][i - 1]], score=score)

    def score(self, letter1, letter2):
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

    def calculate_similarity_matrix(self, seq1, seq2):
        """
        >>> nw = NeedlemanWunsch()
        >>> nw.calculate_similarity_matrix("AATC","AACT")
        >>> nw.scoring_matrix
        array([[ 0., -1., -2., -3., -4.],
               [-1.,  4.,  3.,  2.,  1.],
               [-2.,  3.,  8.,  7.,  6.],
               [-3.,  2.,  7.,  7., 12.],
               [-4.,  1.,  6., 16., 15.]])

        Function which calculates the scoring matrix using needleman-wunsch.
        :param seq1: First sequence.
        :param seq2: Second sequence
        :return: void
        """
        seq1 = "€" + seq1
        seq2 = "€" + seq2
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
                if self.scoring_type == ScoringType.DISTANCE:
                    optimum = min(match, delete, insert)
                else:
                    optimum = max(match, delete, insert)
                self.scoring_matrix[i][j] = optimum
                # We then update the traceback matrix.
                cell = None
                if optimum == match:
                    cell = self.create_traceback_cell(cell=cell, operation=Operation.MATCH, traceback_i=i - 1,
                                                      traceback_j=j - 1, score_i=i, score_j=j)
                if optimum == delete:
                    cell = self.create_traceback_cell(cell=cell, operation=Operation.DELETION, traceback_i=i - 1,
                                                      traceback_j=j, score_i=i, score_j=j)

                if optimum == insert:
                    cell = self.create_traceback_cell(cell=cell, operation=Operation.INSERTION, traceback_i=i,
                                                      traceback_j=j - 1, score_i=i, score_j=j)

                self.traceback_matrix[i][j] = cell

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
        traceback_id = 0
        self.explore(self.desired_traceback, traceback_id)

    def generate_alignments(self, sequence1, sequence2):
        """
        Function which creates all optimal alignments.
        :param sequence1:
        :param sequence2:
        :return:
        """
        alignments = []
        for id, traceback in self.tracebacks.items():
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
            alignments.append(
                    Alignment(sequence1="".join(alignment[0][1:]), sequence2="".join(alignment[1][1:]),
                              operations=alignment[2][1:],
                              score=self.desired_traceback.score))
        self.alignments = alignments
        print(self.alignments)

    def run(self, fasta_files, complete_traceback=False):
        """
        >>> blossum62 = MatrixInfo.blosum62
        >>> nw = NeedlemanWunsch(substitution_matrix=blossum62,similarity=True)
        >>> fasta = ["../data/test1.fn", "../data/test2.fn"]
        >>> nw.run(fasta_files=fasta)
        ('test1', Seq('MNSERSDVTLYQPFLDYAIAYMR', SingleLetterAlphabet()), 'test2', Seq('MNSERSDVTLY', SingleLetterAlphabet()), 36.0, [(NSERSDVTLYQPFLDYAIAYMR, NSERSDVTL-----Y-------, 36), (NSERSDVTLYQPFLD, NSERSDVTLY-----, 36)])

        :param fasta_files:
        :param complete_traceback:
        :return:
        """
        sequences = parse_fasta_files(fasta_files)
        seq1 = sequences[0]
        seq2 = sequences[1]
        self.calculate_similarity_matrix(seq1.seq, seq2.seq)
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
    nw = NeedlemanWunsch(match_scoring=1, indel_scoring=-1, mismatch_scoring=-1, substitution_matrix=blossum62,
                         similarity=True)
    fasta_files = ["../data/test1.fn", "../data/test2.fn"]
    res = nw.run(fasta_files=fasta_files)
    print(res)
