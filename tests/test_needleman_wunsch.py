import pprint
import unittest

import numpy as np
from Bio.SubsMat import MatrixInfo

from needleman_wunsch import NeedlemanWunsch


class NeedlemanWunschTest(unittest.TestCase):
    def test_init(self):
        nw = NeedlemanWunsch()
        nw.init_scoring_matrix("AAAC", "AAAC")
        assert np.array_equal(nw.scoring_matrix, np.array([[0., -6., -12., -18., -24.],
                                                           [-6., 0., 0., 0., 0.],
                                                           [-12., 0., 0., 0., 0.],
                                                           [-18., 0., 0., 0., 0.],
                                                           [-24., 0., 0., 0., 0.]]))

    def test_scoring_blossum(self):
        """Testing score calculation using Blossum62 + Gap Penalty = 6"""

        print("######### Testing calculation of scoring matrix. ###########")
        nw = NeedlemanWunsch(substitution_matrix=MatrixInfo.blosum62, gap_penalty=6)
        print("############# Case 1 ##############")
        print("############# START ##############")
        seq1 = "A"
        seq2 = "A"
        nw.calculate_scoring_matrix(seq1, seq2)
        print("SEQ1: %s" % seq1)
        print("SEQ2: %s" % seq2)
        print("RESULT:\n %s" % nw.scoring_matrix)
        np.testing.assert_array_equal(nw.scoring_matrix, np.array([[0., - 6.],
                                                                   [-6., 4.]]))
        print("############# FINISH ##############")
        print("############# Case 2 ##############")
        print("############# START ##############")
        seq1 = "A"
        seq2 = "AT"
        nw.calculate_scoring_matrix(seq1, seq2)
        print("SEQ1: %s" % seq1)
        print("SEQ2: %s" % seq2)
        print("RESULT:\n %s" % nw.scoring_matrix)
        np.testing.assert_array_equal(nw.scoring_matrix, np.array([[0., - 6., -12.],
                                                                   [-6., 4., -2.]]))
        print("############# FINISH ##############")
        print("############# Case 3 ##############")
        print("############# START ##############")
        seq1 = "ILDMDVVEGSAARFDCKVEGYPDPEVMWFKDDNPVKESRHFQIDYDEEGN"
        seq2 = "RDPVKTHEGWGVMLPCNPPAHYPGLSYRWLLNEFPNFIPTDGRHFVSQTT"
        nw.calculate_scoring_matrix(seq1, seq2)
        print("SEQ1: %s" % seq1)
        print("SEQ2: %s" % seq2)
        print("RESULT:\n %s" % pprint.pformat(nw.scoring_matrix))
        self.assertAlmostEqual(nw.scoring_matrix[-1][-1], 4.0)
        print("############# FINISH ##############")

    def test_scoring_pam(self):
        """Testing score calculation using PAM250 + Gap Penalty = 8"""

        print("######### Testing calculation of scoring matrix. ###########")
        nw = NeedlemanWunsch(substitution_matrix=MatrixInfo.pam250, gap_penalty=8)
        print("############# Case 1 ##############")
        print("############# START ##############")
        seq1 = "A"
        seq2 = "A"
        nw.calculate_scoring_matrix(seq1, seq2)
        print("SEQ1: %s" % seq1)
        print("SEQ2: %s" % seq2)
        print("RESULT:\n %s" % nw.scoring_matrix)
        np.testing.assert_array_equal(nw.scoring_matrix, np.array([[0., - 8.],
                                                                   [-8., 2.]]))
        print("############# FINISH ##############")
        print("############# Case 2 ##############")
        print("############# START ##############")
        seq1 = "A"
        seq2 = "AT"
        nw.calculate_scoring_matrix(seq1, seq2)
        print("SEQ1: %s" % seq1)
        print("SEQ2: %s" % seq2)
        print("RESULT:\n %s" % nw.scoring_matrix)
        np.testing.assert_array_equal(nw.scoring_matrix, np.array([[0., - 8., -16.],
                                                                   [-8., 2., -6.]]))
        print("############# FINISH ##############")
        print("############# Case 3 ##############")
        print("############# START ##############")
        seq1 = "ILDMDVVEGSAARFDCKVEGYPDPEVMWFKDDNPVKESRHFQIDYDEEGN"
        seq2 = "RDPVKTHEGWGVMLPCNPPAHYPGLSYRWLLNEFPNFIPTDGRHFVSQTT"
        nw.calculate_scoring_matrix(seq1, seq2)
        print("SEQ1: %s" % seq1)
        print("SEQ2: %s" % seq2)
        print("RESULT:\n %s" % pprint.pformat(nw.scoring_matrix))
        self.assertAlmostEqual(nw.scoring_matrix[-1][-1], 31.0)
        print("############# FINISH ##############")
