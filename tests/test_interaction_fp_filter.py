import unittest
import sys
from rdkit import Chem
import numpy as np

# import scripts
sys.path.insert(1, '/home/sabsr3/xchem/fragment_network_merges/scripts')
from interaction_fp_filter import InteractionFPFilter

mol = 'tests/x0107-x0678t-0.minimised.mol'
fragmentA = 'tests/Mpro-x0107_0A.mol'
fragmentB = 'tests/Mpro-x0678_0A.mol'
protein = 'tests/Mpro-x0107_0A_apo-desolv.pdb'

class TestInteractionFPFilter(unittest.TestCase):

    def test_make_fp(self):
        test_case = InteractionFPFilter(mol, fragmentA, fragmentB, protein)
        test_mol = test_case.get_mol(mol)
        lst = [1, 0, 0, 1, 0, 0, 0, 0]
        self.assertListEqual(list(test_case.make_fp(test_mol)), lst)
    
    def test_similarity_filter(self):
        test_case = InteractionFPFilter(mol, fragmentA, fragmentB, protein)
        result = 'pass'
        self.assertEqual(test_case.similarity_filter(), result)

if __name__ == '__main__':
    unittest.main()
