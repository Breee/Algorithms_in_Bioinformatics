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

import pprint
import unittest

import numpy as np
from Bio.SubsMat import MatrixInfo

from needleman_wunsch import NeedlemanWunsch, ScoringSettings
from utility.utils import parse_fasta_files


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
        nw = NeedlemanWunsch(ScoringSettings(substitution_matrix=MatrixInfo.blosum62, gap_penalty=6))
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
        nw = NeedlemanWunsch(ScoringSettings(substitution_matrix=MatrixInfo.pam250, gap_penalty=8))
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

    def test_run_pam(self):
        pairs_to_result = {('ILDMDVVEGSAARFDCKVEGYPDPEVMWFKDDNPVKESRHFQIDYDEEGN',
                            'RDPVKTHEGWGVMLPCNPPAHYPGLSYRWLLNEFPNFIPTDGRHFVSQTT'): 31,
                           ('ILDMDVVEGSAARFDCKVEGYPDPEVMWFKDDNPVKESRHFQIDYDEEGN',
                            'ISDTEADIGSNLRWGCAAAGKPRPMVRWLRNGEPLASQNRVEVLA')     : 44,
                           ('ILDMDVVEGSAARFDCKVEGYPDPEVMWFKDDNPVKESRHFQIDYDEEGN',
                            'RRLIPAARGGEISILCQPRAAPKATILWSKGTEILGNSTRVTVTSD')    : 13,
                           ('RDPVKTHEGWGVMLPCNPPAHYPGLSYRWLLNEFPNFIPTDGRHFVSQTT',
                            'ISDTEADIGSNLRWGCAAAGKPRPMVRWLRNGEPLASQNRVEVLA')     : 15,
                           ('RDPVKTHEGWGVMLPCNPPAHYPGLSYRWLLNEFPNFIPTDGRHFVSQTT',
                            'RRLIPAARGGEISILCQPRAAPKATILWSKGTEILGNSTRVTVTSD')    : 16,
                           ('ISDTEADIGSNLRWGCAAAGKPRPMVRWLRNGEPLASQNRVEVLA',
                            'RRLIPAARGGEISILCQPRAAPKATILWSKGTEILGNSTRVTVTSD')    : 45}

        sequence_file = '../data/guideline_tests/needlemanwunsch.fa'
        sequences = parse_fasta_files([sequence_file])
        # init the needleman
        settings = ScoringSettings(substitution_matrix=MatrixInfo.pam250, gap_penalty=8,
                                   similarity=True)
        nw = NeedlemanWunsch(settings, complete_traceback=False, verbose=False)
        results = nw.pairwise_alignments(sequences)
        for result in results:
            seqs = (str(result.seq1.seq), str(result.seq2.seq))
            expected_score = pairs_to_result[seqs]
            self.assertEqual(result.score, expected_score)

    def test_run_blossum(self):
        pairs_to_result = {('ILDMDVVEGSAARFDCKVEGYPDPEVMWFKDDNPVKESRHFQIDYDEEGN',
                            'RDPVKTHEGWGVMLPCNPPAHYPGLSYRWLLNEFPNFIPTDGRHFVSQTT'): 4,
                           ('ILDMDVVEGSAARFDCKVEGYPDPEVMWFKDDNPVKESRHFQIDYDEEGN',
                            'ISDTEADIGSNLRWGCAAAGKPRPMVRWLRNGEPLASQNRVEVLA')     : 37,
                           ('ILDMDVVEGSAARFDCKVEGYPDPEVMWFKDDNPVKESRHFQIDYDEEGN',
                            'RRLIPAARGGEISILCQPRAAPKATILWSKGTEILGNSTRVTVTSD')    : -4,
                           ('RDPVKTHEGWGVMLPCNPPAHYPGLSYRWLLNEFPNFIPTDGRHFVSQTT',
                            'ISDTEADIGSNLRWGCAAAGKPRPMVRWLRNGEPLASQNRVEVLA')     : 3,
                           ('RDPVKTHEGWGVMLPCNPPAHYPGLSYRWLLNEFPNFIPTDGRHFVSQTT',
                            'RRLIPAARGGEISILCQPRAAPKATILWSKGTEILGNSTRVTVTSD')    : 9,
                           ('ISDTEADIGSNLRWGCAAAGKPRPMVRWLRNGEPLASQNRVEVLA',
                            'RRLIPAARGGEISILCQPRAAPKATILWSKGTEILGNSTRVTVTSD')    : 24}

        sequence_file = '../data/guideline_tests/needlemanwunsch.fa'
        sequences = parse_fasta_files([sequence_file])
        # init the needleman
        settings = ScoringSettings(substitution_matrix=MatrixInfo.blosum62, gap_penalty=6,
                                   similarity=True)
        nw = NeedlemanWunsch(settings, complete_traceback=False, verbose=False)
        results = nw.pairwise_alignments(sequences)
        for result in results:
            seqs = (str(result.seq1.seq), str(result.seq2.seq))
            expected_score = pairs_to_result[seqs]
            self.assertEqual(result.score, expected_score)


if __name__ == '__main__':
    suite = unittest.TestSuite()
    suite.addTest(NeedlemanWunschTest())
    unittest.TextTestRunner().run(suite)
