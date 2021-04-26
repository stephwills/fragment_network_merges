"""Tests the find merges script"""
import unittest
from scripts.find_merges import *

# update/check tests to reflect change in script to filter synthons
# test cases - two fragments to merge and a fake smiles
names = ['test_x0107', 'test_x1382', 'fail_case']
smiles = ['CC(=O)Nc1cnccc1C', 'CC(NC(=O)CCl)c1cccc(Cl)c1', 'CC(NC(=O)CCl)c1cccc(Cl)c1s']

class TestNode(unittest.TestCase):
    """Tests the node-related functions"""

    def test_filter_for_nodes(self):
        """Tests the smiles are filtered for those in the network"""
        test_case = filter_for_nodes(smiles, names)
        smiles_filtered = ['CC(=O)Nc1cnccc1C', 'CC(NC(=O)CCl)c1cccc(Cl)c1']
        names_filtered = ['test_x0107', 'test_x1382']
        self.assertEqual(test_case[0], smiles_filtered)
        self.assertEqual(test_case[1], names_filtered)

    def test_get_combinations(self):
        """Tests the combinations of fragments are generated correctly"""
        # test case
        filtered_smiles, filtered_names = filter_for_nodes(smiles, names)
        # actual answer
        fragment_pairs = [('CC(=O)Nc1cnccc1C', 'CC(NC(=O)CCl)c1cccc(Cl)c1'),
                          ('CC(NC(=O)CCl)c1cccc(Cl)c1', 'CC(=O)Nc1cnccc1C')]
        name_pairs = [('test_x0107', 'test_x1382'),
                      ('test_x1382', 'test_x0107')]
        test_case = get_combinations(filtered_smiles, filtered_names)
        self.assertListEqual(test_case[0], fragment_pairs)
        self.assertListEqual(test_case[1], name_pairs)

class TestMerge(unittest.TestCase):
    """Tests the merge-related functions"""

    def test_get_expansions(self):
        """Tests the function generates synthons and merges correctly"""
        # create test case
        filtered_smiles, filtered_names = filter_for_nodes(smiles, names)
        fragment_pairs, name_pairs = get_combinations(filtered_smiles, filtered_names)
        pair1_smiles, pair1_names = fragment_pairs[0], name_pairs[0]
        # run expansion
        test_results = get_expansions(pair1_smiles, pair1_names)
        # get the number of synthons and smiles generated
        num_synthons = len(test_results)
        num_smiles = 0
        for synthon in test_results:
            num_smiles += len(test_results[synthon])
        self.assertEqual(num_synthons, 7)
        self.assertEqual(num_smiles, 3751)

if __name__ == '__main__':
    unittest.main()
