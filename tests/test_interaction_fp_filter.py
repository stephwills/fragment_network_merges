"""Tests the interaction fp script"""
import unittest
import numpy as np
from rdkit import Chem
from scripts.interaction_fp_filter import InteractionFPFilter

# create test cases
mol = 'tests/x0107-x0678t-0.minimised.mol'
fragmentA = 'tests/Mpro-x0107_0A.mol'
fragmentB = 'tests/Mpro-x0678_0A.mol'
protein = 'tests/Mpro-x0107_0A_apo-desolv.pdb'

class TestInteractionFPFilter(unittest.TestCase):
    """Tests the InteractionFPFilter class"""

    def test_make_fp(self):
        """Tests the interaction fp is correctly calculated"""
        test_case = InteractionFPFilter(mol, fragmentA, fragmentB, protein)
        test_mol = test_case.get_mol(mol)
        lst = [1, 0, 0, 1, 0, 0, 0, 0]
        self.assertListEqual(list(test_case.make_fp(test_mol)), lst)
    
    def test_similarity_filter(self):
        """Checks that molecules correctly pass and fail the filter"""
        test_case = InteractionFPFilter(mol, fragmentA, fragmentB, protein)
        result = 'pass'
        self.assertEqual(test_case.similarity_filter(), result)

if __name__ == '__main__':
    unittest.main()
