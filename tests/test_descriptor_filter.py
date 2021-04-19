"""Tests the descriptor filter script"""
import unittest
from scripts.descriptor_filter import calculate_properties, descriptor_filter

class TestDescriptorFilter(unittest.TestCase):
    """Tests the DescriptorFilter class"""

    def test_calculate_properties(self):
        """Checks that a molecule's properties and number of violations are calculated correctly"""
        violations = calculate_properties('NC(=O)CN1CCCOc2cccc(c2)CCC2CCCCN2CC1')
        self.assertEqual(violations, 1)

    def test_filter(self):
        """Checks that molecules correctly pass and fail the filter"""
        passing_smiles = 'NC(=O)CN1CCCOc2cccc(c2)CCC2CCCCN2CC1'
        failing_smiles = 'CC(C)C1=NN=C(C)N1[C@H]1C[C@@H]2CC[C@H](C1)N2CC[C@H](NC(=O)C1CCC(F)(F)CC1)C1=CC=CC=C1'
        passing_case = descriptor_filter(passing_smiles)
        failing_case = descriptor_filter(failing_smiles)
        self.assertEqual(passing_case, 'pass')
        self.assertEqual(failing_case, 'fail')

if __name__ == '__main__':
    unittest.main()
