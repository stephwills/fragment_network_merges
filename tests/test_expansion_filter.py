"""Tests the expansion filter script"""

import unittest

from filter.expansion_filter import ExpansionFilter
from merge.preprocessing import get_mol

fragmentA_path = get_mol("nsp13", "x0276_0B", False)
fragmentB_path = get_mol("nsp13", "x0034_0B", False)
fragmentA = get_mol("nsp13", "x0276_0B", True)
fragmentB = get_mol("nsp13", "x0034_0B", True)


class TestExpansionFilter(unittest.TestCase):
    """Tests the expansion filter function"""

    def test_expansion_filter_failing(self):
        """
        Tests the filter fails where most of the synthon is also in fragment A
        (so looks like an expansion).
        """
        smi = "CC(C)NCC(C)(C)CN(C)C(=O)c1ccccc1F"
        synthon = "Fc1ccccc1[Xe]"
        filter = ExpansionFilter([smi], [synthon], fragmentA_path, fragmentB_path)
        failing_case = filter.filter_smi(smi, synthon, fragmentA, fragmentB)
        self.assertEqual(failing_case, False)

    def test_expansion_filter_passing(self):
        """Tests the filter passes molecules correctly"""
        smi = 'Fc1ccc(CCNc2ccc(F)cc2)cc1'
        synthon = "Fc1ccccc1[Xe]"
        filter = ExpansionFilter([smi], [synthon], fragmentA_path, fragmentB_path)
        failing_case = filter.filter_smi(smi, synthon, fragmentA, fragmentB)
        self.assertEqual(failing_case, True)


if __name__ == "__main__":
    unittest.main()
