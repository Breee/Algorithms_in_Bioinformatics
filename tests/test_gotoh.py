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

import math
import pprint
import unittest

import numpy as np
from Bio.SubsMat import MatrixInfo

from gotoh import Gotoh
from utility.utils import parse_fasta_files


class GotohTest(unittest.TestCase):
    def test_init(self):
        got = Gotoh()
        got.init_scoring_matrices("AAAC", "AAAC")
        print(got.scoring_matrix_D)
        print(got.scoring_matrix_P)
        print(got.scoring_matrix_Q)
        assert np.array_equal(got.scoring_matrix_D, np.array([[0., -12., -13., -14., -15.],
                                                              [-12., 0., 0., 0., 0.],
                                                              [-13., 0., 0., 0., 0.],
                                                              [-14., 0., 0., 0., 0.],
                                                              [-15., 0., 0., 0., 0.]]))
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
        """Testing score calculation using Blossum62 + Gap Penalty = -11, gap extend = -1"""

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

    def test_run_pam(self):

        pairs_to_result = {('ILDMDVVEGSAARFDCKVEGYPDPEVMWFKDDNPVKESRHFQIDYDEEGN',
                            'RDPVKTHEGWGVMLPCNPPAHYPGLSYRWLLNEFPNFIPTDGRHFVSQTT'): 33,
                           ('ILDMDVVEGSAARFDCKVEGYPDPEVMWFKDDNPVKESRHFQIDYDEEGN',
                            'ISDTEADIGSNLRWGCAAAGKPRPMVRWLRNGEPLASQNRVEVLA'):      60,
                           ('ILDMDVVEGSAARFDCKVEGYPDPEVMWFKDDNPVKESRHFQIDYDEEGN',
                            'RRLIPAARGGEISILCQPRAAPKATILWSKGTEILGNSTRVTVTSD'):     30,
                           ('RDPVKTHEGWGVMLPCNPPAHYPGLSYRWLLNEFPNFIPTDGRHFVSQTT',
                            'RRLIPAARGGEISILCQPRAAPKATILWSKGTEILGNSTRVTVTSD'):     9,
                           ('ISDTEADIGSNLRWGCAAAGKPRPMVRWLRNGEPLASQNRVEVLA',
                            'RRLIPAARGGEISILCQPRAAPKATILWSKGTEILGNSTRVTVTSD'):     41,

                           ('RDPVKTHEGWGVMLPCNPPAHYPGLSYRWLLNEFPNFIPTDGRHFVSQTT',
                            'ISDTEADIGSNLRWGCAAAGKPRPMVRWLRNGEPLASQNRVEVLA'):      17
                           }
        sequence_file = '../data/guideline_tests/needlemanwunsch.fa'
        sequences = parse_fasta_files([sequence_file])
        gt = Gotoh(substitution_matrix=MatrixInfo.pam250,
                   gap_penalty=11,
                   gap_extend=1,
                   similarity=True,
                   verbose=False, complete_traceback=True)
        results = gt.pairwise_alignments(sequences)
        for result in results:
            seqs = (str(result.seq1), str(result.seq2))
            expected_score = pairs_to_result[seqs]
            self.assertEqual(result.score, expected_score)
            # print(len(result.alignments))
        # import numpy
        # with numpy.printoptions(threshold=numpy.inf):
        #    print(gt.scoring_matrix_D)

    def test_run_blossum(self):
        pairs_to_result = {('ILDMDVVEGSAARFDCKVEGYPDPEVMWFKDDNPVKESRHFQIDYDEEGN',
                            'RDPVKTHEGWGVMLPCNPPAHYPGLSYRWLLNEFPNFIPTDGRHFVSQTT'): 0,
                           ('ILDMDVVEGSAARFDCKVEGYPDPEVMWFKDDNPVKESRHFQIDYDEEGN',
                            'ISDTEADIGSNLRWGCAAAGKPRPMVRWLRNGEPLASQNRVEVLA'):      41,
                           ('ILDMDVVEGSAARFDCKVEGYPDPEVMWFKDDNPVKESRHFQIDYDEEGN',
                            'RRLIPAARGGEISILCQPRAAPKATILWSKGTEILGNSTRVTVTSD'):     5,
                           ('RDPVKTHEGWGVMLPCNPPAHYPGLSYRWLLNEFPNFIPTDGRHFVSQTT',
                            'RRLIPAARGGEISILCQPRAAPKATILWSKGTEILGNSTRVTVTSD'):     -4,
                           ('ISDTEADIGSNLRWGCAAAGKPRPMVRWLRNGEPLASQNRVEVLA',
                            'RRLIPAARGGEISILCQPRAAPKATILWSKGTEILGNSTRVTVTSD'):     18,

                           ('RDPVKTHEGWGVMLPCNPPAHYPGLSYRWLLNEFPNFIPTDGRHFVSQTT',
                            'ISDTEADIGSNLRWGCAAAGKPRPMVRWLRNGEPLASQNRVEVLA'):      -5
                           }

        sequence_file = '../data/guideline_tests/needlemanwunsch.fa'
        sequences = parse_fasta_files([sequence_file])
        gt = Gotoh(substitution_matrix=MatrixInfo.blosum62,
                   gap_penalty=11,
                   gap_extend=1,
                   similarity=True,
                   verbose=False, complete_traceback=True)
        results = gt.pairwise_alignments(sequences)
        for result in results:
            seqs = (str(result.seq1), str(result.seq2))
            expected_score = pairs_to_result[seqs]
            self.assertEqual(result.score, expected_score)
            print(len(result.alignments))


if __name__ == '__main__':
    suite = unittest.TestSuite()
    suite.addTest(GotohTest())
    unittest.TextTestRunner().run(suite)
