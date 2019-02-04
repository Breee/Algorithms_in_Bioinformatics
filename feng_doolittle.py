#!usr/bin/python3
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

from logger.log import setup_custom_logger
from utility.utils import MultiAlignment, Result, count_gaps_in_pairwise_alignment, \
    count_occurences_symbol_in_word, parse_input, replace_with_neutral_symbol

LOGGER = setup_custom_logger("feng", logfile="feng.log")

import argparse
import math
from needleman_wunsch import NeedlemanWunsch
from xpgma import Xpgma, Node, GuideTree
import copy
import difflib
from Bio.SubsMat import MatrixInfo
from utility.utils import Clustering


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

    def __init__(self, match_scoring=1, indel_scoring=-1, mismatch_scoring=-1, gap_penalty=6,
                 substitution_matrix=MatrixInfo.blosum62, verbose=False, clustering_method=Clustering.UPGMA):
        self.match_scoring = match_scoring
        self.indel_scoring = indel_scoring
        self.mismatch_scoring = mismatch_scoring
        # Subsitution matrices (PAM/BLOSSOM), default is blossom62.
        self.substitution_matrix = substitution_matrix
        self.verbose = verbose
        self.clustering_method = clustering_method

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

        s_eff = (s_ab.score - s_rand) / (s_max - s_rand)

        # negative values make no sense.
        if s_eff <= 0.0:
            score = 1
        else:
            score = - math.log(s_eff)
        LOGGER.info("New score: %.5f" % score)
        return score

    def compute_best_alignment_one_to_many(self, leaf: Node, alignment: MultiAlignment):
        """
        Function which finds the best alignment, by calculating alignments between a sequence and many sequences.
        :param leaf: Node object which is a leaf
        :param alignment: MultiAlignment object
        :return: alignment, index of best alignment, alignment score.
        """
        assert leaf.is_leaf()
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

    def compute_best_alignment_many_to_many(self, alignment1: MultiAlignment, alignment2: MultiAlignment):
        """
        Function which finds the best alignment, by calculating alignment between two lists of sequences.
        :param alignment1: MultiAlignment object
        :param alignment2: MultiAlignment object
        :return: best_alignment, index in alignment1, index in alignment2, best_score, overall_score
        """
        best_alignment = None
        index1 = None
        index2 = None
        best_score = None
        overall_score = 0
        sequences1 = alignment1.sequences
        sequences2 = alignment2.sequences
        nw = NeedlemanWunsch()
        for i, seq1 in enumerate(sequences1):
            for j, seq2 in enumerate(sequences2):
                result = nw.run(seq1, seq2)
                if best_score is None or result.score > best_score:
                    best_score = result.score
                    best_alignment = result.alignments[0]
                    index1 = i
                    index2 = j
                # the score is the addition of all pairwise scores.
                overall_score += result.score
        return [best_alignment.sequence1, best_alignment.sequence2], index1, index2, best_score, overall_score

    def operation1(self, leaf1: Node, leaf2: Node) -> MultiAlignment:
        """
        Compute best pairwise alignment,
        change occurences of gap symbol to X
        :param leaf1: Node object
        :param leaf2: Node object
        :return: MultiAlignment object

        >>> from Bio.SeqRecord import SeqRecord
        >>> feng = FengDoolittle()
        >>> res = feng.operation1(leaf1=Node(sequence=SeqRecord("AAACGA"),name=None, cost=None),\
                            leaf2=Node(sequence=SeqRecord("AAA"), name=None,cost=None))
        >>> res.sequences[0].seq
        'AAACGA'
        >>> res.sequences[1].seq
        'XAAXXA'
        """
        assert leaf1.is_leaf() and leaf2.is_leaf()
        nw = NeedlemanWunsch()
        result = nw.run(leaf1.sequence, leaf2.sequence)
        multi_alignment = MultiAlignment(sequences=[result.alignments[0].sequence1, result.alignments[0].sequence2],
                                         score=result.score)
        multi_alignment.sequences = replace_with_neutral_symbol(multi_alignment.sequences)
        return multi_alignment

    def operation2(self, leaf: Node, alignment: MultiAlignment) -> MultiAlignment:
        """
        Align a sequence S to an alignment A
        1. Compute pairwise alignment score between S and all sequences S' in A
        2. Align S according to the sequence in the best pairwise alignment
        3. Change occurences of gap symbol to X
        :param leaf:
        :param alignment:
        :return:

        >>> from Bio.SeqRecord import SeqRecord
        >>> feng = FengDoolittle()
        >>> leaf = Node(None,None,sequence=SeqRecord("AA"))
        >>> alignment = MultiAlignment(sequences=[SeqRecord("AAXX"),SeqRecord("GGCC")], score=4)
        >>> res = feng.operation2(leaf,alignment)
        >>> res.sequences[0].seq
        'AAXX'
        >>> res.sequences[1].seq
        'AAXX'
        >>> res.sequences[2].seq
        'GGCC'
        >>> res.score
        0.0
        """
        best_alignment, index, score = self.compute_best_alignment_one_to_many(leaf, alignment)
        best_alignment = replace_with_neutral_symbol(best_alignment)
        alignment.sequences.pop(index)
        alignment.sequences.insert(index, best_alignment[0])
        alignment.sequences.insert(index, best_alignment[1])
        alignment.score += score
        return alignment

    @staticmethod
    def reforge_with_gaps(sequences, blueprint_sequence_index, old):
        """
        Function which will update sequences and inserts gaps at the positions, where the new one has them.
        :param sequences:
        :param blueprint_sequence_index:
        :return:
        >>> from Bio.SeqRecord import SeqRecord
        >>> sequences = [SeqRecord("AA"), SeqRecord("XAXA"), SeqRecord("XXGG")]
        >>> res = FengDoolittle.reforge_with_gaps(sequences, blueprint_sequence_index=1, old=SeqRecord("AA"))
        >>> res[0].seq
        'XAXA'
        >>> res[1].seq
        'XAXA'
        >>> res[2].seq
        'XXGG'
        """
        new = sequences[blueprint_sequence_index]
        add_list = []
        # first we have to find out where new gaps have been inserted.
        # if old and new are equal size, we assume nothing changed.
        if len(new) == len(old):
            return sequences
        else:
            for i, s in enumerate(difflib.ndiff(old, new)):
                if s[0] == ' ':
                    continue
                elif s[0] == '+':
                    # print(u'Add "{}" to position {}'.format(s[-1], i))
                    add_list.append(i)
        for seq in sequences:
            if seq.seq != new.seq:
                split = list(seq.seq)
                for el in add_list:
                    if split[el] != "X":
                        split.insert(el, "X")
                seq.seq = "".join(split)
        return sequences

    def operation3(self, alignment1: MultiAlignment, alignment2: MultiAlignment) -> MultiAlignment:
        """
        Align an alignment A1 to an alignment A2.
        For each pair of sequences S1 in A1 and S2 in A2 compute pairwise alignment score
        Align A1 and A2 according to the pairwise alignment with minimal distance.
        Change occurences of gap symbol to X

        :param alignment1: MultiAlignment object
        :param alignment2: MultiAlignment object
        :return: MultiAlignment object
        """
        best_alignment, index1, index2, score, overall_score = self.compute_best_alignment_many_to_many(alignment1,
                                                                                                        alignment2)
        best_alignment = replace_with_neutral_symbol(best_alignment)
        # remove old elements.
        alignment1.sequences.pop(index1)
        alignment1.sequences.insert(index1, best_alignment[0])
        old = alignment2.sequences.pop(index2)
        alignment2.sequences.insert(index2, best_alignment[1])
        alignment2.sequences = self.reforge_with_gaps(alignment2.sequences, index2, old)
        # merge sequences lists
        new_sequences = alignment1.sequences + alignment2.sequences
        # add new best_alignment.
        new_alignment = MultiAlignment(sequences=new_sequences, score=overall_score)
        return new_alignment

    def traverse(self, node: Node) -> MultiAlignment:
        """
        Core function of the algorithm, it traverses a guide tree node to the bottom
        and returns a MultiAlignment object.
        There are three cases which are handled:
        Case 1: Apply operation 1 on the sequences returned by the leafs and return the alignment
        Case 2: Apply operation 2 on the sequence and the alignment and return the alignment
        Case 3: Both Children are inner nodes, Apply operation 3 on the alignments and return the alignment
        :param node: Node object
        :return: MultiAlignment object
        """
        if node.is_leaf():
            return node.sequence
        else:
            child1 = node.children[0]
            child2 = node.children[1]
            # Case 1: Both children are leaf nodes
            # Apply operation 1 on the sequences returned by the leafs and return the alignment
            if child1.is_leaf() and child2.is_leaf():
                return self.operation1(child1, child2)
            # Case 2: One Child is leaf and one node is inner node
            # Apply operation 2 on the sequence and the alignment and return the
            # alignment
            elif child1.is_leaf():
                return self.operation2(child1, self.traverse(child2))
            elif child2.is_leaf():
                return self.operation2(child2, self.traverse(child1))
            # Case 3: Both Children are inner nodes
            # Apply operation 3 on the alignments and return the alignment
            else:
                return self.operation3(self.traverse(child1), self.traverse(child2))

    def compute_msa(self, guidetree: GuideTree) -> MultiAlignment:
        """
        Function to kickoff msa calculation
        :param guidetree:
        :return: MultiAlignment object
        """
        msa = self.traverse(guidetree.root)
        return msa

    def run(self, sequences):
        """
        Run function for feng doolittle.
        :param sequences: a list of SeqRecords
        :return: MultiAlignment object
        """
        # perform pairwise sequence alignments
        nw = NeedlemanWunsch(verbose=self.verbose)
        alignments = nw.pairwise_alignments(sequences)
        LOGGER.info("Needleman Wunsch Alignments:\n%s" % "\n".join([str(x) for x in alignments]))
        # Convert the scores to approximate pairwise evolutionary distances.
        for alignment in alignments:
            alignment.score = self.convert_to_evolutionary_distances(alignment)
        # 2. Construct a guide tree
        # init the xpgma
        xpgma = Xpgma(clustering_method=self.clustering_method)
        tree = xpgma.run(alignments)
        # 3. Start from the root of the tree to compute MSA.
        msa = self.compute_msa(tree)
        res_str = "\n".join([x.seq for x in msa.sequences])
        LOGGER.info(f'Tree: {tree}')
        LOGGER.info("GENERATED MSA:\nSCORE:%f\nMSA:\n\n%s" % (msa.score, res_str))
        return msa


def process_program_arguments():
    LOGGER.info(f'Arguments:\n {args}')
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
    if args.clustering_mode.lower() == 'upgma':
        args.clustering_mode = Clustering.UPGMA
    elif args.clustering_mode.lower() == 'wpgma':
        args.clustering_mode = Clustering.WPGMA
    else:
        LOGGER.critical(
                "UNKNOWN parameter for clustering mode. Choose UPGMA or WPGMA, falling back to UPGMA")
        args.clustering_mode = Clustering.UPGMA


def run_feng_doolittle():
    sequences = parse_input(args.input, args.file_filter)
    if len(sequences) < 2:
        LOGGER.warn("We received not enough sequences. Make sure you called the program correctly.")
        exit(1)
    elif len(sequences) >= 2:
        feng = FengDoolittle(substitution_matrix=args.substitution_matrix, gap_penalty=args.gap_penalty,
                             verbose=args.verbose, clustering_method=args.clustering_mode)
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
    parser.add_argument('-m', '--match', type=float, default=1.0,
                        help='match score')
    parser.add_argument('-mm', '--mismatch', type=float, default=-1.0,
                        help='mismatch score')
    parser.add_argument('-g', '--gap_penalty', type=float, default=6.0,
                        help='gap penalty')
    parser.add_argument('-s', '--substitution_matrix', type=str, default="BLOSSUM",
                        help='Substitution Matrix (BLOSSUM | PAM | NONE) default is BLOSSUM')
    parser.add_argument('-cm', '--clustering_mode', type=str, default="UPGMA",
                        help='UPGMA | WPGMA')

    args = parser.parse_args()
    main()
else:
    pass
