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

from unittest import TestCase

from needleman_wunsch import NeedlemanWunsch
from utility import utils
from xpgma import Xpgma


class TestXpgma(TestCase):
    def test_create_distance_matrix(self):
        nw = NeedlemanWunsch()
        sequences = utils.parse_fasta_files(["../data/xpgma/xpgma1.fa"])
        alignments = nw.pairwise_alignments(sequences)
        xpgma = Xpgma()
        xpgma.create_distance_matrix(alignments)
        self.assertDictEqual(xpgma.distances,
                             {'A': {'B': 4.0, 'C': 4.0}, 'B': {'A': 4.0, 'C': 4.0}, 'C': {'A': 4.0, 'B': 4.0}})

    def test_get_minimum_distance(self):
        min = Xpgma.get_minimum_distance({'A': {'B': 4.0, 'C': 0.0}, 'B': {'C': -3.0}})
        expected = (-3.0, 'B', 'C')
        self.assertEqual(min, expected)

    def test_get_dist(self):
        xpgma = Xpgma()
        xpgma.distances = {'A': {'B': 4.0, 'C': 0.0}, 'B': {'C': -3.0}}
        d1 = xpgma.get_dist("A", "C")
        e1 = 0.0
        self.assertEqual(d1, e1)

        d2 = xpgma.get_dist("B", "C")
        e2 = -3.0
        self.assertEqual(d2, e2)

        d3 = xpgma.get_dist("A", "B")
        e3 = 4.0
        self.assertEqual(d3, e3)

        d4 = xpgma.get_dist("B", "A")
        e4 = 4.0
        self.assertEqual(d4, e4)

        d5 = xpgma.get_dist("C", "A")
        e5 = 0.0
        self.assertEqual(d5, e5)

        d6 = xpgma.get_dist("C", "B")
        e6 = -3.0
        self.assertEqual(d6, e6)

        d7 = xpgma.get_dist("A", "A")
        e7 = 0.0
        self.assertEqual(d7, e7)

    def test_merge_clusters(self):
        xpgma = Xpgma()
        xpgma.distances = {'A': {'B': 4.0, 'C': 4.0}, 'B': {'A': 4.0, 'C': 4.0}, 'C': {'A': 4.0, 'B': 4.0}}
        xpgma.merge_clusters('A', 'B')
        actual = xpgma.distances
        expected = {'AB': {'C': 4.0}, 'C': {'AB': 4.0}}
        self.assertDictEqual(actual, expected)

    def test_calculate_guide_tree(self):
        nw = NeedlemanWunsch()
        sequences = utils.parse_fasta_files(["../data/xpgma/xpgma1.fa"])
        alignments = nw.pairwise_alignments(sequences)
        xpgma = Xpgma()
        xpgma.create_distance_matrix(alignments)
        guidetree = xpgma.calculate_guide_tree()
        expected = '((A:2.00,B:2.00):0.00,C:2.00)'
        expected_nodes = "{'A': A:2.00, 'B': B:2.00, 'C': C:2.00, 'AB': (A:2.00,B:2.00):0.00, 'ABC': ABC:NONE}"
        self.assertEqual(str(guidetree), expected)
        self.assertEqual(str(guidetree.nodes), expected_nodes)

    def test_add_to_guidetree(self):
        self.fail()
