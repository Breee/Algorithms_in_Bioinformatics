from unittest import TestCase

from feng_doolittle import FengDoolittle
from needleman_wunsch import NeedlemanWunsch
from utility import utils

TESTFILES = [""]


class TestFengDoolittle(TestCase):
    def test_convert_to_evolutionary_distances(self):
        # perform pairwise sequence alignments
        nw = NeedlemanWunsch()
        sequences = utils.parse_fasta_files(["../data/feng_test/conversion.fa"])
        alignments = nw.pairwise_alignments(sequences)
        feng = FengDoolittle()
        # Convert the scores to approximate pairwise evolutionary distances.
        alignment = alignments[0]
        print(f'Alignment: {alignment} ')
        alignment.score = feng.convert_to_evolutionary_distances(alignment)
        print(f'Score: {alignment.score} ')
        self.assertAlmostEqual(first=2.70805020110221, second=alignment.score)

    def test_run(self):
        # perform pairwise sequence alignments
        sequences = utils.parse_fasta_files(["../data/feng_test/feng1.fa"])
        feng = FengDoolittle()
        res = feng.run(sequences)
        expected_order = ["ILDMDVVEGSAARFDCKVEGYPDPEVMWFKDDNPVKESRHFQIDYDEEGN",
                          "ISDTEADIGSNLRWGCXAAAGKPRPMVRWLRNGEPLXASQNXRVEVXXLAX",
                          "RRLIPAARGGEISILCQPRAAXPKATILWSKXGTEILGNXSTRVTVXTXSD",
                          "RDPVKTHEGWGVMLPCNPPAHYPGLSYRWLLNEFPNFIPTDGRHFVSQXTT"]
        print(res.sequences)
        print(res.score)
