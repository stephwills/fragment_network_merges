"""Tests the preprocessing script"""

import unittest
import numpy as np
from scripts.preprocessing import check_fragment_pairs, get_distance, get_distance_between_fragments, get_mol, get_smiles

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

    def test_get_distance(self):
        """Tests distance calculated correctly"""
        point_A = np.array([0,0,0])
        point_B = np.array([0,8,6])
        dist = 10.0
        self.assertEqual(dist, get_distance(point_A, point_B))

    def test_get_distance_between_fragments(self):
        """Tests the distance between fragments calculated OK"""
        test_fA = get_mol('Mpro', 'x0107_0A')
        test_fB = get_mol('Mpro', 'x0434_0A')
        dist = 0.39528976713292285
        self.assertEqual(dist, get_distance_between_fragments(test_fA, test_fB))
    
    def test_check_fragment_pairs(self):
        """Tests that fragments that are far apart are filtered out"""
        # 2nd pair is fail case
        fragments = [('x11427_0A', 'x10870_0A'),  ('x2659_0A', 'x10876_0A')]
        smiles = [(get_smiles('Mpro', f[0]), get_smiles('Mpro', f[1])) for f in fragments]
        self.assertEqual(len(check_fragment_pairs(smiles, fragments, 'Mpro')), 2)

if __name__ == '__main__':
    unittest.main()
