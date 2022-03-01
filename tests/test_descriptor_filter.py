"""Tests the descriptorvfilter script"""

import unittest

from filter.descriptor_filter import DescriptorFilter


class TestDescriptorFilter(unittest.TestCase):
    """Tests the DescriptorFilter class"""

    def test_calculate_properties(self):
        """Checks that a molecule's properties and number of violations are calculated correctly"""
        smi = "CC(C)C1=NN=C(C)N1[C@H]1C[C@@H]2CC[C@H](C1)N2CC[C@H](NC(=O)C1CCC(F)(F)CC1)C1=CC=CC=C1"
        filter = DescriptorFilter([smi])
        violations = filter.calculate_properties(smi)
        self.assertEqual(violations, 2)

    def test_filter(self):
        """Checks that molecules correctly pass and fail the filter"""
        smis = [
            "NC(=O)CN1CCCOc2cccc(c2)CCC2CCCCN2CC1",
            "CC(C)C1=NN=C(C)N1[C@H]1C[C@@H]2CC[C@H](C1)N2CC[C@H](NC(=O)C1CCC(F)(F)CC1)C1=CC=CC=C1",
        ]
        filter = DescriptorFilter(smis)
        actual_results = [True, False]
        results = filter.filter_all()[0]
        self.assertEqual(results, actual_results)


if __name__ == "__main__":
    unittest.main()
