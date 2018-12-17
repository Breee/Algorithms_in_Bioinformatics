import math
import pprint
import unittest

import numpy as np
from Bio.SubsMat import MatrixInfo

from gotoh import Gotoh


class GotohTest(unittest.TestCase):
    def test_init(self):
        got = Gotoh()
        got.init_scoring_matrices("AAAC", "AAAC")
        print(got.scoring_matrix_D)
        print(got.scoring_matrix_P)
        print(got.scoring_matrix_Q)
        assert np.array_equal(got.scoring_matrix_D, np.array([[0., -9., -12., -15., -18.],
                                                              [-9., 0., 0., 0., 0.],
                                                              [-12., 0., 0., 0., 0.],
                                                              [-15., 0., 0., 0., 0.],
                                                              [-18., 0., 0., 0., 0.]]))
        assert np.array_equal(got.scoring_matrix_P, np.array([[-math.inf, -math.inf, -math.inf, -math.inf, -math.inf],
                                                              [-math.inf, 0., 0., 0., 0.],
                                                              [-math.inf, 0., 0., 0., 0.],
                                                              [-math.inf, 0., 0., 0., 0.],
                                                              [-math.inf, 0., 0., 0., 0.]]))
        assert np.array_equal(got.scoring_matrix_Q, np.array([[-math.inf, -math.inf, -math.inf, -math.inf, -math.inf],
                                                              [-math.inf, 0., 0., 0., 0.],
                                                              [-math.inf, 0., 0., 0., 0.],
                                                              [-math.inf, 0., 0., 0., 0.],
                                                              [-math.inf, 0., 0., 0., 0.]]))

    def test_scoring_blossum(self):
        """Testing score calculation using Blossum62 + Gap Penalty = 6"""

        print("######### Testing calculation of scoring matrix. ###########")
        got = Gotoh(substitution_matrix=MatrixInfo.blosum62, gap_penalty=-11, gap_extend=-1)
        print("############# Case 1 ##############")
        print("############# START ##############")
        seq1 = "A"
        seq2 = "A"
        got.calculate_scoring_matrix(seq1, seq2)
        print("SEQ1: %s" % seq1)
        print("SEQ2: %s" % seq2)
        print("RESULT:\n %s" % got.scoring_matrix_D)
        np.testing.assert_array_equal(got.scoring_matrix_D, np.array([[0., -12.0],
                                                                      [-12., 4.]]))
        print("############# FINISH ##############")
        print("############# Case 2 ##############")
        print("############# START ##############")
        seq1 = "A"
        seq2 = "AT"
        got.calculate_scoring_matrix(seq1, seq2)
        print("SEQ1: %s" % seq1)
        print("SEQ2: %s" % seq2)
        print("RESULT:\n %s" % got.scoring_matrix_D)
        np.testing.assert_array_equal(got.scoring_matrix_D, np.array([[0., -12., -13.],
                                                                      [-12., 4., -8.]]))
        print("############# FINISH ##############")
        print("############# Case 3 ##############")
        print("############# START ##############")
        seq1 = "RDPVKTHEGWGVMLPCNPPAHYPGLSYRWLLNEFPNFIPTDGRHFVSQTT"
        seq2 = "ISDTEADIGSNLRWGCAAAGKPRPMVRWLRNGEPLASQNRVEVLA"
        got.calculate_scoring_matrix(seq1, seq2)
        print("SEQ1: %s" % seq1)
        print("SEQ2: %s" % seq2)
        print("RESULT:\n %s" % pprint.pformat(got.scoring_matrix_D))
        self.assertAlmostEqual(got.scoring_matrix_D[-1][-1], -5.0)
        print("############# FINISH ##############")

    def test_scoring_pam(self):
        """Testing score calculation using PAM250 + Gap Penalty = -11, gap extend = -1"""

        print("######### Testing calculation of scoring matrix. ###########")
        got = Gotoh(substitution_matrix=MatrixInfo.pam250, gap_penalty=-11, gap_extend=-1)
        print("############# Case 1 ##############")
        print("############# START ##############")
        seq1 = "A"
        seq2 = "A"
        got.calculate_scoring_matrix(seq1, seq2)
        print("SEQ1: %s" % seq1)
        print("SEQ2: %s" % seq2)
        print("RESULT:\n %s" % got.scoring_matrix_D)
        np.testing.assert_array_equal(got.scoring_matrix_D, np.array([[0., -12.],
                                                                      [-12., 2.]]))
        print("############# FINISH ##############")
        print("############# Case 2 ##############")
        print("############# START ##############")
        seq1 = "A"
        seq2 = "AT"
        got.calculate_scoring_matrix(seq1, seq2)
        print("SEQ1: %s" % seq1)
        print("SEQ2: %s" % seq2)
        print("RESULT:\n %s" % got.scoring_matrix_D)
        np.testing.assert_array_equal(got.scoring_matrix_D, np.array([[0., - 12., -13.],
                                                                      [-12., 2., -10.]]))
        print("############# FINISH ##############")
        print("############# Case 3 ##############")
        print("############# START ##############")
        seq1 = "ILDMDVVEGSAARFDCKVEGYPDPEVMWFKDDNPVKESRHFQIDYDEEGN"
        seq2 = "RDPVKTHEGWGVMLPCNPPAHYPGLSYRWLLNEFPNFIPTDGRHFVSQTT"
        got.calculate_scoring_matrix(seq1, seq2)
        print("SEQ1: %s" % seq1)
        print("SEQ2: %s" % seq2)
        print("RESULT:\n %s" % pprint.pformat(got.scoring_matrix_D))
        self.assertAlmostEqual(got.scoring_matrix_D[-1][-1], 33)
        print("############# FINISH ##############")


if __name__ == '__main__':
    suite = unittest.TestSuite()
    suite.addTest(GotohTest())
    unittest.TextTestRunner().run(suite)
