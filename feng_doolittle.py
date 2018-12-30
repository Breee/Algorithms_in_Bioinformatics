from logger.log import setup_custom_logger
from utility.utils import Result, count_gaps_in_pairwise_alignment, count_occurences_symbol_in_word, parse_input

LOGGER = setup_custom_logger("feng", logfile="feng.log")

import argparse
import math
from needleman_wunsch import NeedlemanWunsch
from xpgma import Xpgma
import copy


class FengDoolittle(object):
    """
    Class which implements the feng doolittle
    1. Calculate the pairwise alignment scores, and convert them to distances.
    2. Use a clustering algorithm to construct a tree from the distances.
    3. Traverse the nodes in their order of addition to the tree, repeatedly align the child nodes (sequences or
    alignments).
    Features of this heuristic:
    - Highest scoring pairwise alignment determines the alignment to two groups.
    - "Once a gap, always a gap": replace gaps in alignments by a neutral character.
    """

    def __init__(self, verbose=False):
        pass

    def convert_to_evolutionary_distances(self, pairwise_alignment_result: Result) -> float:
        """Converts similarity score from a pairwise alignment to a distance score
        using approximation algorithm

        D(a,b) = - log(S_{a,b}^{eff})
        S_{a,b}^{eff} = (S(a,b) - S_{rand}) / (S_{a,b}^{max} - S_{rand})
        S_{rand} = (1/|A|) * (sum_{x,y in \Sigma \times \Sigma} S(x,y) * N_a(x) * N_b(y)) + gaps(A) * S(-,*)
        S_{a,b}^{max} = (S(a,a) + S(b,b)) / 2
        """
        alignment = pairwise_alignment_result.alignments[0]
        LOGGER.info("Converting similarity to evolutionary distances.")
        LOGGER.info("Alignment: %s" % alignment)

        seq1 = copy.deepcopy(alignment.sequence1)
        seq1.seq = seq1.seq.replace("-", "")

        seq2 = copy.deepcopy(alignment.sequence2)
        seq2.seq = seq2.seq.replace("-", "")

        nw = NeedlemanWunsch()
        s_ab = nw.run(seq1, seq2)
        s_aa = nw.run(seq1, seq1)
        s_bb = nw.run(seq2, seq2)

        s_max = (s_aa.score + s_bb.score) / 2
        s_rand = (1 / len(alignment.sequence1)) * \
                 sum([nw.score(nw.alphabet.letters[i],
                               nw.alphabet.letters[j])
                      * count_occurences_symbol_in_word(seq1.seq, nw.alphabet.letters[i])
                      * count_occurences_symbol_in_word(seq2.seq, nw.alphabet.letters[j])
                      for i in range(len(nw.alphabet.letters)) for j in range(len(nw.alphabet.letters))]) \
                 + count_gaps_in_pairwise_alignment(alignment) * nw.gap_penalty

        # prevent division by zero.
        if s_max == s_rand:
            s_rand = s_rand - 0.0001

        S_eff = (s_ab.score - s_rand) / (s_max - s_rand)

        # negative values make no sense.
        if S_eff <= 0.0:
            score = 1
        else:
            score = - math.log(S_eff)
        LOGGER.info("New score: %.5f" % score)
        return score

    def compute_best_alignment_one_to_many(self, leaf: Node, alignment: MultiAlignment):
        best_alignment = None
        index = None
        best_score = None
        leaf_sequence = leaf.sequence
        sequences = alignment.sequences
        nw = NeedlemanWunsch()
        for i, seq in enumerate(sequences):
            result = nw.run(leaf_sequence, seq)
            if best_score is None or result.score > best_score:
                best_score = result.score
                best_alignment = result.alignments[0]
                index = i
        return [best_alignment.sequence1, best_alignment.sequence2], index, best_score
    def run(self, sequences):
        # init the xpgma
        # perform pairwise sequence alignments
        nw = NeedlemanWunsch(verbose=args.verbose)
        alignments = nw.pairwise_alignments(sequences)
        LOGGER.info("Needleman Wunsch Alignments:\n%s" % "\n".join([str(x) for x in alignments]))
        # Convert the scores to approximate pairwise evolutionary distances.
        for alignment in alignments:
            alignment.score = self.convert_to_evolutionary_distances(alignment)
        # 2. Construct a guide tree
        xpgma = Xpgma()
        tree = xpgma.run(alignments)
        # 3. Start from the first node that has been added to the guide tree and align the child nodes
        
        # 4. Repeat step 3. For all other nodes in the order in which they were added to the tree.
        # Do this until all sequences have been aligned.


def process_program_arguments():
    if not args.input:
        LOGGER.critical(
                "Error, provide input file/files/directory/directories using -i / --input. -h/--help for all "
                "arguments.")
        exit(1)


def run_feng_doolittle():
    sequences = parse_input(args.input, args.file_filter)
    if len(sequences) < 2:
        LOGGER.warn("We received not enough sequences. Make sure you called the program correctly.")
        exit(1)
    elif len(sequences) >= 2:
        feng = FengDoolittle()
        feng.run(sequences)


def main():
    process_program_arguments()
    run_feng_doolittle()
    exit(1)


if __name__ == '__main__':
    # command line arguments.
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input', type=str, required=True, default='',
                        nargs='+',
                        help='A file, a directory or multiple directories. directories are processed recursively.')
    parser.add_argument('--file-filter', type=str, default='',
                        help='A regex to define a filter, which will be applied to the files parsed.')
    parser.add_argument('-v', '--verbose', action='store_true', default=False,
                        help='verbose output.')

    args = parser.parse_args()
    main()
else:
    pass
