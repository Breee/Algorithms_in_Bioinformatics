from needleman_wunsch import NeedlemanWunsch
from utility.utils import parse_fasta_files


def test_seq_1_2():
    """Example testing the dummy implementation."""

    nw = NeedlemanWunsch()
    fasta_files = ["../data/test3.fa"]
    sequences = parse_fasta_files(fasta_files)
    print(sequences)
    seq1 = sequences[0]
    seq2 = sequences[1]
    result = nw.run(seq1, seq2, complete_traceback=False)
    (id_seq1, seq1, id_seq2, seq2, score, alignments) = result
    print((id_seq1, seq1, id_seq2, seq2, score, alignments))
