"""Tests the expansion filter script"""

import os
import unittest

from filter.expansion_filter import ExpansionFilter
from merge.preprocessing import get_mol
from rdkit import Chem

frag_dir = os.path.join("tests", "test_Fragalysis")
fragmentA_path = get_mol("nsp13", "x0176_0B", False, frag_dir)
fragmentB_path = get_mol("nsp13", "x0034_0B", False, frag_dir)
fragmentA = get_mol("nsp13", "x0176_0B", True, frag_dir)
fragmentB = get_mol("nsp13", "x0034_0B", True, frag_dir)


class TestExpansionFilter(unittest.TestCase):
    """Tests the expansion filter function"""

    def test_expansion_filter_volume(self):
        """
        Tests the filter fails where most of the synthon is also in fragment A
        (so looks like an expansion).
        """
        smi = "CCc1ccccc1CCN[S](C)(=O)=O"
        synthon = "Fc1ccccc1[Xe]"
        filter = ExpansionFilter([smi], [synthon], fragmentA_path, fragmentB_path)
        failing_case = filter.filter_smi(smi, fragmentA, fragmentB, 0.9)
        self.assertEqual(failing_case, True)


    def test_expansion_filter_failing(self):
        """
        Tests the filter fails where most of the synthon is also in fragment A
        (so looks like an expansion).
        """
        smi = "CS(=O)(=O)NCCC1=CC=CC=C1Cl"
        synthon = "Fc1ccccc1[Xe]"
        filter = ExpansionFilter([smi], [synthon], fragmentA_path, fragmentB_path)
        failing_case = filter.filter_smi(smi, fragmentA, fragmentB, 0.9)
        self.assertEqual(failing_case, False)


    def test_expansion_filter_passing(self):
        """Tests the filter passes molecules correctly"""
        smi = 'CS(=O)(=O)NCCC1=CC=CC=C1F'
        synthon = "Fc1ccccc1[Xe]"
        filter = ExpansionFilter([smi], [synthon], fragmentA_path, fragmentB_path)
        failing_case = filter.filter_smi(smi, fragmentA, fragmentB, 0.9)
        self.assertEqual(failing_case, True)


if __name__ == "__main__":
    unittest.main()
