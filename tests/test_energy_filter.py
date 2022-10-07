"""Tests the descriptor filter script"""

import os
import unittest
from unittest.mock import patch
from rdkit import Chem

from filter.energy_filter import EnergyFilter, main

test_pass_sdf = os.path.join("tests", "test_data", "energy_filter_pass_mols.sdf")
test_fail_sdf = os.path.join("tests", "test_data", "energy_filter_fail_mols.sdf")
output_sdf = os.path.join("tests", "test_data", "test_energy_filter_output.sdf")


class TestEnergyFilter(unittest.TestCase):
    """Tests the DescriptorFilter class"""

    def test_pass_mols(self):
        """Checks pass mols"""
        mols = [i for i in Chem.SDMolSupplier(test_pass_sdf)]
        filter = EnergyFilter(mols=mols)
        res = filter.filter_smi(mols[0], n_conf=50, energy_threshold=7)
        self.assertEqual(res, True)

    def test_fail_mols(self):
        """Checks fail mols"""
        mols = [i for i in Chem.SDMolSupplier(test_fail_sdf)]
        filter = EnergyFilter(mols=mols)
        res = filter.filter_smi(mols[0], n_conf=50, energy_threshold=7)
        self.assertEqual(res, False)

    def test_main(self):
        """Check that main executes correctly and produces output file"""
        with patch(
            "sys.argv",
            [
                "filter/energy_filter.py",
                "-i",
                test_fail_sdf,
                "-o",
                output_sdf,
                "-e",
                "11",
                "-n",
                "55",
            ],
        ):
            main()
        self.assertTrue(os.path.exists(output_sdf))
        os.remove(output_sdf)


if __name__ == "__main__":
    unittest.main()
