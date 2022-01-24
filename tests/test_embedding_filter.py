"""Tests the embedding filter script"""
import os
import unittest
import numpy as np

from rdkit import Chem

from filter.embedding_filter import EmbeddingFilter, remove_xe
from merge.preprocessing import get_mol


frag_dir = os.path.join('tests', 'test_Fragalysis')
fragmentA_path = get_mol('Mpro', 'x0107_0A', frag_dir)
fragmentB_path = get_mol('Mpro', 'x0678_0A', frag_dir)
fragmentA = get_mol('Mpro', 'x0107_0A', frag_dir, True)
fragmentB = get_mol('Mpro', 'x0678_0A', frag_dir, True)
passing_smi = 'NC(=O)CN1CCC2(C1)CC1(C2)OCCO1'
failing_smi = 'CC(C)(C)C(=O)C=C1SCC(=O)N1CC(N)=O'
synthon = 'NC(=O)C[Xe]'
passing_mol = Chem.MolFromSmiles(passing_smi)
failing_mol = Chem.MolFromSmiles(failing_smi)
synthon_mol = Chem.MolFromSmiles(synthon)
filter = EmbeddingFilter([passing_smi, failing_smi], [synthon, synthon], fragmentA_path, fragmentB_path)

print(filter.filter_all()[0])

class TestEmbeddingFilter(unittest.TestCase):
    """Tests the embedding filter functions"""

    def test_get_mcs(self):
        """Checks that the MCS is correctly identified"""
        mcs = filter._get_mcs(passing_mol, fragmentA)
        self.assertEqual(Chem.MolToSmarts(mcs), '[#6](-[#6])-[#7]-[#6]:,-[#6](:,-[#6]:,-[#6])-[#6]')

    def test_get_distance(self):
        """Checks the distance is calculated correctly"""
        actual_distance = (np.sqrt(((1-4)**2) + ((2-5)**2) + ((3-6)**2)))
        test_distance = filter.get_distance(np.array([1,2,3]), np.array([4,5,6]))
        self.assertEqual(actual_distance, test_distance)

    def test_remove_xe(self):
        """Checks the xenon is removed correctly"""
        xe_removed = 'CC(N)=O'
        self.assertEqual(Chem.MolToSmiles(remove_xe(synthon_mol)), xe_removed)

    def test_calc_energy(self):
        """Tests the energy calculated correctly"""
        energy = 81.32052954057191
        test_unconstrained_energy = filter.calc_unconstrained_energy(passing_mol)
        self.assertEqual(energy, test_unconstrained_energy)

    def test_filter_passing(self):
        """Checks that the filter correctly identifies molecules that can be embedded"""
        passing_case = filter.filter_smi(passing_smi, fragmentA, fragmentB, synthon)
        self.assertEqual(passing_case[0], True)

    def test_filter_failing(self):
        """Checks that the filter correctly identifies molecules that can be embedded"""
        failing_case = filter.filter_smi(failing_smi, fragmentA, fragmentB, synthon)
        self.assertEqual(failing_case[0], False)


if __name__ == '__main__':
    unittest.main()
