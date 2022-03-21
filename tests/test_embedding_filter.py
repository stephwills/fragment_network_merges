"""Tests the embedding filter script"""
import os
import unittest
import numpy as np

from rdkit import Chem
from filter.embedding_filter import EmbeddingFilter, remove_xe
from merge.preprocessing import get_mol


fragmentA_path = get_mol("nsp13", "x0212_0B")
fragmentB_path_passing = get_mol("nsp13", "x0438_0B")
fragmentB_path_failing = get_mol("nsp13", "x0311_0B")
fragmentA = get_mol("nsp13", "x0212_0B", True)
fragmentB_passing = get_mol("nsp13", "x0438_0B", True)
fragmentB_failing = get_mol("nsp13", "x0311_0B", True)
smi = "CN(CCCC(=O)Nc1cnn(-c2ccccc2)c1)S(=O)(=O)c1ccc(F)cc1"
synthon = "[Xe]c1ccccc1"
mol = Chem.MolFromSmiles(smi)
synthon_mol = Chem.MolFromSmiles(synthon)
filter = EmbeddingFilter([smi], [synthon], fragmentA_path, fragmentB_path_passing)


class TestEmbeddingFilter(unittest.TestCase):
    """Tests the embedding filter functions"""

    def test_get_mcs(self):
        """Checks that the MCS is correctly identified"""
        mol = Chem.MolFromSmiles("NC(=O)CN1CCC2(C1)CC1(C2)OCCO1")
        mcs = filter._get_mcs(mol, fragmentA)
        self.assertEqual(Chem.MolToSmarts(mcs), "[#6](-&!@[#6])-&!@[#7]")

    def test_get_distance(self):
        """Checks the distance is calculated correctly"""
        actual_distance = np.sqrt(((1 - 4) ** 2) + ((2 - 5) ** 2) + ((3 - 6) ** 2))
        test_distance = filter.get_distance(np.array([1, 2, 3]), np.array([4, 5, 6]))
        self.assertEqual(actual_distance, test_distance)

    def test_remove_xe(self):
        """Checks the xenon is removed correctly"""
        synthon = "NC(=O)C[Xe]"
        synthon_mol = Chem.MolFromSmiles(synthon)
        xe_removed = "CC(N)=O"
        self.assertEqual(Chem.MolToSmiles(remove_xe(synthon_mol)), xe_removed)

    def test_calc_energy(self):
        """Tests the energy calculated correctly"""
        energy = 72.43
        test_unconstrained_energy = round(filter.calc_unconstrained_energy(mol), 2)
        self.assertEqual(energy, test_unconstrained_energy)

    def test_filter_failing(self):
        """Checks that the filter correctly identifies molecules that can be embedded"""
        filter = EmbeddingFilter(
            [smi], [synthon], fragmentA_path, fragmentB_path_passing
        )
        passing_case = filter.filter_smi(smi, fragmentA, fragmentB_passing, synthon)
        self.assertEqual(passing_case[0], False)

    def test_filter_passing(self):
        """Checks that the filter correctly identifies molecules that can be embedded"""
        filter = EmbeddingFilter(
            [smi], [synthon], fragmentA_path, fragmentB_path_failing
        )
        failing_case = filter.filter_smi(smi, fragmentA, fragmentB_failing, synthon)
        self.assertEqual(failing_case[0], True)


if __name__ == "__main__":
    unittest.main()
