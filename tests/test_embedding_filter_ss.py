"""Tests the embedding filter script"""

import unittest

from filter.embedding_filter_ss import EmbeddingFilterSS
from utils.utils import get_mol

fragmentA = get_mol("nsp13", "x0212_0B", True)
fragmentB_passing = get_mol("nsp13", "x0438_0B", True)
fragmentB_failing = get_mol("nsp13", "x0311_0B", True)
smi = "CN(CCCC(=O)Nc1cnn(-c2ccccc2)c1)S(=O)(=O)c1ccc(F)cc1"


class TestEmbeddingFilter(unittest.TestCase):
    """Tests the embedding filter functions"""

    def test_filter_failing(self):
        """Checks that the filter correctly identifies molecules that can be embedded"""
        filter = EmbeddingFilterSS()
        passing_case = filter.filter_smi(
            smi, fragmentA, fragmentB_passing, 10.0, 50, 0.5
        )
        self.assertEqual(passing_case[0], False)

    def test_filter_passing(self):
        """Checks that the filter correctly identifies molecules that can be embedded"""
        filter = EmbeddingFilterSS()
        failing_case = filter.filter_smi(
            smi, fragmentA, fragmentB_failing, 10.0, 50, 0.5
        )
        self.assertEqual(failing_case[0], False)


if __name__ == "__main__":
    unittest.main()
