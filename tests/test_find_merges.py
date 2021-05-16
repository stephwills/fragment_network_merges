"""Tests the find merges script"""
import unittest
from rdkit import Chem
from rdkit.Chem import rdmolfiles
from scripts.find_merges import filter_for_nodes, filter_synthons, get_combinations, get_expansions, get_synthons, substructure_check
from scripts.preprocessing import get_mol

class TestFindMerges(unittest.TestCase):
    """Tests the node-related functions"""

    def test_filter_for_nodes(self):
        """Tests the smiles are filtered for those in the network"""
        # create test case
        names = ['test_x0107', 'test_x1382', 'fail_case']
        smiles = ['CC(=O)Nc1cnccc1C', 'CC(NC(=O)CCl)c1cccc(Cl)c1', 'CC(NC(=O)CCl)c1cccc(Cl)c1s']
        test_case = filter_for_nodes(smiles, names)
        # filtered answer
        smiles_filtered = ['CC(=O)Nc1cnccc1C', 'CC(NC(=O)CCl)c1cccc(Cl)c1']
        names_filtered = ['test_x0107', 'test_x1382']
        self.assertEqual(test_case[0], smiles_filtered)
        self.assertEqual(test_case[1], names_filtered)

    def test_get_combinations(self):
        """Tests the combinations of fragments are generated correctly"""
        # create test case
        names = ['test_x0107', 'test_x1382']
        smiles = ['CC(=O)Nc1cnccc1C', 'CC(NC(=O)CCl)c1cccc(Cl)c1']
        # actual answer
        fragment_pairs = [('CC(=O)Nc1cnccc1C', 'CC(NC(=O)CCl)c1cccc(Cl)c1'),
                          ('CC(NC(=O)CCl)c1cccc(Cl)c1', 'CC(=O)Nc1cnccc1C')]
        name_pairs = [('test_x0107', 'test_x1382'),
                      ('test_x1382', 'test_x0107')]
        test_case = get_combinations(smiles, names)
        self.assertListEqual(test_case[0], fragment_pairs)
        self.assertListEqual(test_case[1], name_pairs)
    
    def test_get_synthons(self):
        """Tests the function returns all the synthons for a smiles"""
        smiles_11812 = 'O=C(c1cncc2ccccc12)N1CCN(c2ccccc2)CC1'  # test case
        synthons = get_synthons(smiles_11812)
        num_synthons = 18  # correct num synthons
        self.assertEqual(len(synthons), num_synthons)

    def test_filter_synthons(self):
        """Tests that synthons that have <3 carbons are filtered out"""
        # create test cases
        passing_synthon = 'O=C(c1cncc2ccccc12)N1CCN(c2ccccc2)CC1'
        failing_synthon = 'O=C[Xe]'
        self.assertEqual(filter_synthons(passing_synthon), passing_synthon)
        self.assertIsNone(filter_synthons(failing_synthon), None)
    
    def test_substructure_check(self):
        fragment_A = rdmolfiles.MolFromMolFile('tests/Mpro-x0107_0A.mol')
        fragment_B = rdmolfiles.MolFromMolFile('tests/Mpro-x0678_0A.mol')
        failing_synthon = '[Xe]c1ccncc1'
        passing_synthon = 'Cc1ccncc1[Xe]'
        self.assertEqual(substructure_check(passing_synthon, fragment_A, fragment_B), passing_synthon)
        self.assertIsNone(substructure_check(failing_synthon, fragment_A, fragment_B))

        

    # def test_get_expansions(self):
    #     """Tests the function generates synthons and merges correctly"""
    #     # create test case
    #     names = ['test_x0107', 'test_x1382']
    #     smiles = ['CC(=O)Nc1cnccc1C', 'CC(NC(=O)CCl)c1cccc(Cl)c1']
    #     fragment_pairs, name_pairs = get_combinations(smiles, names)
    #     pair1_smiles, pair1_names = fragment_pairs[0], name_pairs[0]
    #     # run expansion
    #     test_results = get_expansions(pair1_smiles, pair1_names)
    #     # get the number of synthons and smiles generated
    #     num_synthons = len(test_results)
    #     num_smiles = 0
    #     for synthon in test_results:
    #         num_smiles += len(test_results[synthon])
    #     self.assertEqual(num_synthons, 7)
    #     self.assertEqual(num_smiles, 3751)

if __name__ == '__main__':
    unittest.main()

# print(len(get_synthons('O=C(c1cncc2ccccc12)N1CCN(c2ccccc2)CC1')))
# print(filter_synthons('[Xe]c1cncc2ccccc12.c1ccccc1'))
# print(filter_synthons('O=C[Xe]'))
# fragment_A = rdmolfiles.MolFromMolFile('tests/Mpro-x0107_0A.mol')
# fragment_B = rdmolfiles.MolFromMolFile('tests/Mpro-x0678_0A.mol')
# smiles_A = Chem.MolToSmiles(fragment_A)
# synthons = get_synthons(smiles_A)
# for s in synthons:
#     print(substructure_check(s, fragment_A, fragment_B))

