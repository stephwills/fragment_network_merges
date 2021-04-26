"""Tests the interaction fp script"""
import unittest
import numpy as np
from rdkit import Chem
from scripts.interaction_fp_filter import *

# create test cases
mol = 'tests/x0107-x0678t-0.minimised.mol'
fragmentA = 'Mpro/aligned/Mpro-x0678_0A/Mpro-x0678_0A.mol'
fragmentB = 'Mpro/aligned/Mpro-x12696_0A/Mpro-x12696_0A.mol'
protein = 'tests/Mpro-x0107_0A_apo-desolv.pdb'

class TestInteractionFPFilter(unittest.TestCase):
    """Tests the interaction fp filter functions"""

    def test_make_fp(self):
        """Tests the interaction fp is correctly calculated"""
        test_mol = get_mol(mol)
        test_protein = get_protein(protein)
        test_fp = make_fp(test_protein, test_mol)
        lst = [1, 0, 0, 0, 0, 0, 0, 0]
        self.assertListEqual(list(test_fp), lst)
    
    def test_similarity_filter(self):
        """Checks that molecules correctly pass and fail the filter"""
        test_case = similarity_filter(mol, fragmentA, fragmentB, protein)
        result = 'fail'
        self.assertEqual(test_case, result)

if __name__ == '__main__':
    unittest.main()
