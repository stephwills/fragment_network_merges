"""Tests the descriptor filter script"""

import os
import unittest
from unittest.mock import patch

from fragment_network_merges.filter import NonringBondFilter, main

test_sdf = os.path.join('tests', 'test_data', 'nonring_bond_filter_mols.sdf')
output_sdf = os.path.join("tests", "test_data", "test_nonring_bond_filter_output.sdf")


class TestDescriptorFilter(unittest.TestCase):
    """Tests the DescriptorFilter class"""

    def test_long_linker(self):
        """Checks a long linker fails the filter"""
        smi = 'C1=CC=CC=C1CCCCCCCCCCCCC2=CC=CC=C2'
        filter = NonringBondFilter([smi])
        res = filter.filter_smi(smi)
        self.assertEqual(False, res)

    def test_long_sidechain(self):
        """Checks that a long sidechain fails the filter"""
        smi = 'C1=CC=CC=C1CCCCCCCCCCCC'
        filter = NonringBondFilter([smi])
        res = filter.filter_smi(smi)
        self.assertEqual(res, False)

    def test_passing(self):
        """Checks a molecule correctly passes the filter"""
        smi = 'C1=CC=CC=C1CCC2=CC=CC=C2'
        filter = NonringBondFilter([smi])
        res = filter.filter_smi(smi)
        self.assertEqual(res, True)

    def test_main(self):
        """Check that main executes correctly and produces output file"""
        with patch(
            "sys.argv",
            [
                "filter/nonring_bond_filter.py",
                "-i",
                test_sdf,
                "-o",
                output_sdf,
                "-s",
                "7"
            ],
        ):
            main()
        self.assertTrue(os.path.exists(output_sdf))
        os.remove(output_sdf)


if __name__ == "__main__":
    unittest.main()
