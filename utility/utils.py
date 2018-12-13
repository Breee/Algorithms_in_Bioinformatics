from enum import Enum

from Bio import SeqIO


class Operation(Enum):
    MATCH = 1,
    INSERTION = 2,
    DELETION = 3,
    MISMATCH = 4,
    ROOT = 0


class ScoringType(Enum):
    SIMILARITY = 0,
    DISTANCE = 1


class Alignment(object):
    def __init__(self, sequence1, sequence2, operations, score):
        self.sequence1 = sequence1
        self.sequence2 = sequence2
        self.operations = operations
        self.score = score
        print(self)

    def __repr__(self):
        return "Alignment: (%s, %s), Score: %d" % (self.sequence1, self.sequence2, self.score)


def parse_fasta_files(files):
    """
    >>> files = ["../data/test1.fn", "../data/test2.fn"]
    >>> parse_fasta_files(files)
    [SeqRecord(seq=Seq('MNSERSDVTLYQPFLDYAIAYMR', SingleLetterAlphabet()), id='test1', name='test1', description=' test1', dbxrefs=[]), SeqRecord(seq=Seq('MNSERSDVTLY', SingleLetterAlphabet()), id='test2', name='test2', description='test2', dbxrefs=[])]

    :param files:
    :return:
    """
    records = []
    for file in files:
        records.extend(list(SeqIO.parse(file, "fasta")))
    return records
