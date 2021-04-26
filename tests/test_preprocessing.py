"""Tests the preprocessing script"""

import unittest
from scripts.preprocessing import *

class TestPreprocessing(unittest.TestCase):
    """Tests the preprocessing functions"""

    def test_get_smiles(self):
        """Tests smiles are retrieved correctly"""
        test_smi = get_smiles('Mpro', 'x0107_0A')
        smi = 'CC(=O)Nc1cnccc1C'
        self.assertEqual(test_smi, smi)
    
    def test_get_mol(self):
        """Tests mol is retrieved correctly"""
        test_mol = get_mol('Mpro', 'x0107_0A')
        # test the 'type' correct
        type_mol = type(get_mol('Mpro', 'x0107_0A')).__name__
        self.assertEqual('Mol', type_mol)
    
    def test_get_distance_between_fragments(self):
        """Tests the distance between fragments calculated OK"""
        test_fA = get_mol('Mpro', 'x0107_0A')
        test_fB = get_mol('Mpro', 'x0434_0A')
        dist = 0.39528976713292285
        self.assertEqual(dist, get_distance_between_fragments(test_fA, test_fB))

if __name__ == '__main__':
    unittest.main()
