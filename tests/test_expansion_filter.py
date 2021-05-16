"""Tests the expansion filter script"""
import unittest
from rdkit import Chem
from scripts.expansion_filter import check_for_mcs, expansion_filter
from scripts.preprocessing import get_mol

class TestExpansionFilter(unittest.TestCase):
    """Tests the expansion filter function"""

    def test_check_for_mcs_passing(self):
        """Tests the MCS is retrieved correctly"""
        fragmentA = get_mol('nsp13', 'x0176_0B')
        fragmentB = get_mol('nsp13', 'x0034_0B')
        merge = Chem.MolFromSmiles('CCC(CCC(O)c1ccccc1F)NC(=O)OC(C)(C)C')
        mols = [fragmentA, fragmentB, merge]
        mcs_smiles = Chem.MolToSmiles(check_for_mcs(mols))
        # MCS is benzene with methyl attached
        correct_mcs = 'Cc1ccccc1'
        self.assertEqual(mcs_smiles, correct_mcs)
    
    def test_check_for_mcs_failing(self):
        """Tests that if no MCS, function returns None"""
        fragmentA = get_mol('nsp13', 'x0176_0B')
        fragmentB = get_mol('nsp13', 'x0034_0B')
        merge = Chem.MolFromSmiles('')
        mols = [fragmentA, fragmentB, merge]
        mcs_mol = check_for_mcs(mols)
        self.assertIsNone(mcs_mol)

    def test_expansion_filter_failing(self):
        """
        Tests the filter fails where most of the synthon is also in fragment A
        (so looks like an expansion).
        """
        fragmentA = get_mol('nsp13', 'x0176_0B')
        fragmentB = get_mol('nsp13', 'x0034_0B')
        merge = Chem.MolFromSmiles('CCC(CCC(O)c1ccccc1F)NC(=O)OC(C)(C)C')
        synthon = Chem.MolFromSmiles('Fc1ccccc1[Xe]')
        failing_case = expansion_filter(merge, fragmentA, fragmentB, synthon)
        self.assertEqual(failing_case, 'fail')
    
    def test_expansion_filter_passing(self):
        """Tests the filter passes molecules correctly"""
        fragmentA = get_mol('nsp13', 'x0176_0B')
        fragmentB = get_mol('nsp13', 'x0438_0B')
        merge = Chem.MolFromSmiles('CCC(=O)NC1CCN(CCNS(C)(=O)=O)C1')
        synthon = Chem.MolFromSmiles('CCC(=O)N[Xe]')
        passing_case = expansion_filter(merge, fragmentA, fragmentB, synthon)
        self.assertEqual(passing_case, 'pass')

if __name__ == '__main__':
    unittest.main()
