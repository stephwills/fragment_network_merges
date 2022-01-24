"""Tests the expansion filter script"""

import os
import unittest
from rdkit import Chem

from filter.expansion_filter import ExpansionFilter
from merge.preprocessing import get_mol

frag_dir = os.path.join('tests', 'test_Fragalysis')
fragmentA_path = get_mol('nsp13', 'x0176_0B', frag_dir)
fragmentB_path = get_mol('nsp13', 'x0034_0B', frag_dir)
fragmentA = get_mol('nsp13', 'x0176_0B', frag_dir, True)
fragmentB = get_mol('nsp13', 'x0034_0B', frag_dir, True)

class TestExpansionFilter(unittest.TestCase):
    """Tests the expansion filter function"""

    def test_check_for_mcs_passing(self):
        """Tests the MCS is retrieved correctly"""
        smi = 'CCC(CCC(O)c1ccccc1F)NC(=O)OC(C)(C)C'
        mol = Chem.MolFromSmiles(smi)
        filter = ExpansionFilter([smi], [], fragmentA_path, fragmentB_path)
        mols = [fragmentA, fragmentB, mol]
        mcs_smiles = Chem.MolToSmiles(filter._check_for_mcs(mols))
        # MCS is benzene with methyl attached
        correct_mcs = 'Cc1ccccc1'
        self.assertEqual(mcs_smiles, correct_mcs)

    def test_check_for_mcs_failing(self):
        """Tests that if no MCS, function returns None"""
        smi = ''
        mol = Chem.MolFromSmiles(smi)
        filter = ExpansionFilter([smi], [], fragmentA_path, fragmentB_path)
        mols = [fragmentA, fragmentB, mol]
        mcs_mol = filter._check_for_mcs(mols)
        self.assertIsNone(mcs_mol)

    def test_expansion_filter_failing(self):
        """
        Tests the filter fails where most of the synthon is also in fragment A
        (so looks like an expansion).
        """
        smi = 'CCC(CCC(O)c1ccccc1F)NC(=O)OC(C)(C)C'
        synthon = 'Fc1ccccc1[Xe]'
        filter = ExpansionFilter([smi], [synthon], fragmentA_path, fragmentB_path)
        failing_case = filter.filter_smi(smi, fragmentA, fragmentB, synthon, 0.9)
        self.assertEqual(failing_case, False)

    def test_expansion_filter_passing(self):
        """Tests the filter passes molecules correctly"""
        smi = 'CCC(=O)NC1CCN(CCNS(C)(=O)=O)C1'
        synthon = 'CCC(=O)N[Xe]'
        filter = ExpansionFilter([smi], [synthon], fragmentA_path, fragmentB_path)
        passing_case = filter.filter_smi(smi, fragmentA, fragmentB, synthon, 0.9)
        self.assertEqual(passing_case, True)

    def test_expansion_filter_edge(self):
        """Test edge case where >1 MCS match"""
        fragmentA_ = get_mol('nsp13', 'x0311_0B', frag_dir, True)
        fragmentB_ = get_mol('nsp13', 'x0438_0B', frag_dir, True)
        fragmentA_path_ = get_mol('nsp13', 'x0311_0B', frag_dir)
        fragmentB_path_ = get_mol('nsp13', 'x0438_0B', frag_dir)
        smi = 'CCC(=O)Nc1ccccc1Cc1ccccc1'
        synthon = 'CCC(=O)N[Xe]'
        filter = ExpansionFilter([smi], [synthon], fragmentA_path_, fragmentB_path_)
        passing_case = filter.filter_smi(smi, fragmentA_, fragmentB_, synthon, 0.9)
        self.assertEqual(passing_case, True)

if __name__ == '__main__':
    unittest.main()
