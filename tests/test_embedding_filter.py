import unittest
import sys
from rdkit import Chem

# import scripts
sys.path.insert(1, '/home/sabsr3/xchem/fragment_network_merges/scripts')
from embedding_filter import EmbeddingFilter

# some test cases
passing_smiles = 'NC(=O)CN1CCC2(C1)CC1(C2)OCCO1'
failing_smiles = 'CC(C)(C)C(=O)C=C1SCC(=O)N1CC(N)=O'
synthon = 'NC(=O)C[Xe]'
fragmentA = Chem.MolFromMolFile('tests/Mpro-x0107_0A.mol')
fragmentB = Chem.MolFromMolFile('tests/Mpro-x0678_0A.mol')

# finish writing unit tests
class TestEmbeddingFilter(unittest.TestCase):

    def test_get_mcs(self):
        passing_case = EmbeddingFilter(passing_smiles, fragmentA, fragmentB, synthon)
        mcs = passing_case.get_mcs(Chem.MolFromSmiles(passing_smiles), fragmentA)
        self.assertEqual(Chem.MolToSmarts(mcs), '[#6](-[#6])-[#7]-[#6]:,-[#6](:,-[#6]:,-[#6])-[#6]')

    def test_filter(self):
        passing_case = EmbeddingFilter(passing_smiles, fragmentA, fragmentB, synthon).filter()
        failing_case = EmbeddingFilter(failing_smiles, fragmentA, fragmentB, synthon).filter()
        self.assertIsNotNone(passing_case)
        self.assertIsNone(failing_case)

if __name__ == '__main__':
    unittest.main()
