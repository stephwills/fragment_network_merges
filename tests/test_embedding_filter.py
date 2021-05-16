"""Tests the embedding filter script"""
import unittest
import numpy as np
from rdkit import Chem
from scripts.embedding_filter import calc_unconstrained_energy, embedding_filter, get_distance, get_mcs, remove_xe

# some test cases
passing_mol = Chem.MolFromSmiles('NC(=O)CN1CCC2(C1)CC1(C2)OCCO1')
failing_mol = Chem.MolFromSmiles('CC(C)(C)C(=O)C=C1SCC(=O)N1CC(N)=O')
synthon = Chem.MolFromSmiles('NC(=O)C[Xe]')
fragmentA = Chem.MolFromMolFile('tests/Mpro-x0107_0A.mol')
fragmentB = Chem.MolFromMolFile('tests/Mpro-x0678_0A.mol')

class TestEmbeddingFilter(unittest.TestCase):
    """Tests the embedding filter functions"""

    def test_get_mcs(self):
        """Checks that the MCS is correctly identified"""
        mcs = get_mcs(passing_mol, fragmentA)
        self.assertEqual(Chem.MolToSmarts(mcs), '[#6](-[#6])-[#7]-[#6]:,-[#6](:,-[#6]:,-[#6])-[#6]')

    def test_get_distance(self):
        """Checks the distance is calculated correctly"""
        actual_distance = (np.sqrt(((1-4)**2) + ((2-5)**2) + ((3-6)**2)))
        test_distance = get_distance(np.array([1,2,3]), np.array([4,5,6]))
        self.assertEqual(actual_distance, test_distance)

    def test_remove_xe(self):
        """Checks the xenon is removed correctly"""
        xe_removed = 'CC(N)=O'
        self.assertEqual(Chem.MolToSmiles(remove_xe(synthon)), xe_removed)

    def test_calc_energy(self):
        """Tests the energy calculated correctly"""
        energy = 81.32052954057191
        test_unconstrained_energy = calc_unconstrained_energy(passing_mol)
        self.assertEqual(energy, test_unconstrained_energy)

    def test_filter(self):
        """Checks that the filter correctly identifies molecules that can be embedded"""
        passing_case = embedding_filter(passing_mol, fragmentA, fragmentB, synthon)
        failing_case = embedding_filter(failing_mol, fragmentA, fragmentB, synthon)
        self.assertEqual(passing_case[1], 'pass')
        self.assertEqual(failing_case[1], 'fail')

if __name__ == '__main__':
    unittest.main()
