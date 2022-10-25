"""
Tests the embedding filter script
"""

import os
import unittest
from unittest.mock import patch

import numpy as np
from filter.embedding_filter import EmbeddingFilter, main
from merge.preprocessing import get_mol
from rdkit import Chem
from utils.filter_utils import calc_unconstrained_energy, get_mcs, remove_xe
from utils.utils import get_distance

fragalysis_dir = "tests/test_Fragalysis"
filter = EmbeddingFilter()
fragmentA = get_mol("nsp13", "x0212_0B", True, fragalysis_dir)
fragmentB_failing = get_mol("nsp13", "x0438_0B", True, fragalysis_dir)
fragmentB_passing = get_mol("nsp13", "x0311_0B", True, fragalysis_dir)
smi = "CN(CCCC(=O)Nc1cnn(-c2ccccc2)c1)S(=O)(=O)c1ccc(F)cc1"
synthon = "[Xe]c1ccccc1"
test_sdf = os.path.join("tests", "test_data", "embedding_filter_mols.sdf")
output_sdf = os.path.join("tests", "test_data", "test_embedding_output.sdf")
test_fragmentA = get_mol("nsp13", "x0034_0B", False, fragalysis_dir)
test_fragmentB = get_mol("nsp13", "x0176_0B", False, fragalysis_dir)


class TestEmbeddingFilter(unittest.TestCase):
    """Tests the embedding filter and associated functions"""

    def test_get_mcs(self):
        """Checks that the MCS is correctly identified"""
        mol = Chem.MolFromSmiles("NC(=O)CN1CCC2(C1)CC1(C2)OCCO1")
        mcs = get_mcs(mol, fragmentA)
        self.assertEqual(Chem.MolToSmarts(mcs), "[#6](-&!@[#6])-&!@[#7]")

    def test_get_distance(self):
        """Checks the distance is calculated correctly"""
        actual_distance = np.sqrt(((1 - 4) ** 2) + ((2 - 5) ** 2) + ((3 - 6) ** 2))
        test_distance = get_distance(np.array([1, 2, 3]), np.array([4, 5, 6]))
        self.assertEqual(actual_distance, test_distance)

    def test_remove_xe(self):
        """Checks the xenon is removed correctly"""
        synthon = "NC(=O)C[Xe]"
        synthon_mol = Chem.MolFromSmiles(synthon)
        xe_removed = "CC(N)=O"
        self.assertEqual(Chem.MolToSmiles(remove_xe(synthon_mol)), xe_removed)

    def test_calc_energy(self):
        """Tests the energy calculated correctly"""
        energy = 71.85
        mol = Chem.MolFromSmiles(smi)
        test_unconstrained_energy = round(calc_unconstrained_energy(mol, 50), 2)
        self.assertEqual(energy, test_unconstrained_energy)

    def test_filter_failing(self):
        """Checks that the filter correctly identifies molecules that can be embedded"""
        failing_case = filter.filter_smi(
            smi, fragmentA, fragmentB_failing, synthon, 7.0, 50, 1.0
        )
        self.assertEqual(failing_case[0], False)

    def test_filter_passing(self):
        """Checks that the filter correctly identifies molecules that can be embedded"""
        passing_case = filter.filter_smi(
            smi, fragmentA, fragmentB_passing, synthon, 7.0, 50, 1.0
        )
        self.assertEqual(passing_case[0], True)

    def test_filter_passing_simsearch(self):
        simsearch_case = filter.filter_smi(
            "O=C(NCCc1ccccc1)c1ccc(O)cc1",
            get_mol("nsp13", "x0208_0A", True, fragalysis_dir),
            get_mol("nsp13", "x0438_0B", True, fragalysis_dir),
            None,
            7.0,
            50,
            1.0,
        )
        self.assertEqual(simsearch_case[0], True)

    def test_main(self):
        """Check that main executes correctly and produces output file"""
        with patch(
            "sys.argv",
            [
                "filter/embedding_filter.py",
                "-i",
                test_sdf,
                "-o",
                output_sdf,
                "-a",
                test_fragmentA,
                "-b",
                test_fragmentB,
                "--n_conformations",
                "55",
            ],
        ):
            main()
        self.assertTrue(os.path.exists(output_sdf))
        os.remove(output_sdf)


if __name__ == "__main__":
    unittest.main()
