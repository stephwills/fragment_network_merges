"""Tests the embedding filter script"""
import unittest
from rdkit import Chem
from scripts.embedding_filter import EmbeddingFilter

# some test cases
passing_smiles = 'NC(=O)CN1CCC2(C1)CC1(C2)OCCO1'
failing_smiles = 'CC(C)(C)C(=O)C=C1SCC(=O)N1CC(N)=O'
synthon = 'NC(=O)C[Xe]'
fragmentA = Chem.MolFromMolFile('tests/Mpro-x0107_0A.mol')
fragmentB = Chem.MolFromMolFile('tests/Mpro-x0678_0A.mol')

# finish writing unit tests
class TestEmbeddingFilter(unittest.TestCase):
    """Tests the EmbeddingFilter class"""

    def test_get_mcs(self):
        """Checks that the MCS is correctly identified"""
        passing_case = EmbeddingFilter(passing_smiles, fragmentA, fragmentB, synthon)
        mcs = passing_case.get_mcs(Chem.MolFromSmiles(passing_smiles), fragmentA)
        self.assertEqual(Chem.MolToSmarts(mcs), '[#6](-[#6])-[#7]-[#6]:,-[#6](:,-[#6]:,-[#6])-[#6]')

    def test_filter(self):
        """Checks that the filter correctly identifies molecules that can be embedded"""
        passing_case = EmbeddingFilter(passing_smiles, fragmentA, fragmentB, synthon).filter()
        failing_case = EmbeddingFilter(failing_smiles, fragmentA, fragmentB, synthon).filter()
        self.assertEqual(passing_case, 'pass')
        self.assertEqual(failing_case, 'fail')

if __name__ == '__main__':
    unittest.main()
