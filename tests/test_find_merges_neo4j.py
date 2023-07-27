"""Tests the find merges neo4j script"""
import os
import unittest

from fragment_network_merges.merge.find_merges_neo4j import MergerFinder_neo4j
from fragment_network_merges.utils import get_mol
from rdkit import Chem

merger = MergerFinder_neo4j()

def cleanup_files(dir):
    files = [os.path.join(dir, f) for f in os.listdir(dir)]
    for file in files:
        os.remove(file)

# results will depend on version of the network -- tests just check there is a result
class TestFindMerges(unittest.TestCase):
    def test_check_for_nodes(self):
        true_smi = Chem.MolToSmiles(Chem.MolFromSmiles("CC=1C=CC(CS(=O)(=O)N)=CC1F"))
        false_smi = "CC(NC(=O)CCl)c1cccc(Cl)c1s"
        true_results = (1, 2)  # number of nodes present, number of fragments queried
        results = merger.check_for_nodes([true_smi, false_smi])
        self.assertEqual(true_results, results)

    def test_get_synthons(self):
        smi = Chem.MolToSmiles(Chem.MolFromSmiles("CC=1C=CC(CS(=O)(=O)N)=CC1F"))
        synthons = merger.get_synthons(smi)
        self.assertTrue(len(synthons) > 0)

    def test_get_expansions(self):
        smi = Chem.MolToSmiles(Chem.MolFromSmiles("CC=1C=CC(CS(=O)(=O)N)=CC1F"))
        smi2 = Chem.MolToSmiles(Chem.MolFromSmiles("CS(=O)(=O)NCCC=1C=CC=CC1"))
        results = merger.get_expansions(
            (smi, smi2), ("x0034_0B", "x0176_0B"), "nsp13", None
        )
        self.assertTrue(len(results) > 0)

    ## new code tests ###
    def test_get_unique_synthons(self):
        nameA = "x0176_0B"
        nameBs = ["x0034_0B", "x0183_0B", "x0438_0B"]
        target = "nsp13"
        unique_synthons, synthon_dict = merger.get_unique_synthons(
            nameA, nameBs, target
        )
        self.assertTrue(len(unique_synthons) > 0)
        self.assertTrue(len(synthon_dict) > 0)
        self.assertTrue(len(synthon_dict["x0034_0B"]) > 0)

    def test_expand_fragmentA(self):
        name_pairs = [["x0034_0B", "x0020_0B"]]
        target = "nsp13"
        output_dir = "tests/test_output"
        merger.expand_fragmentA(name_pairs, target, output_dir=output_dir, working_dir=output_dir, fragalysis_dir='tests/test_Fragalysis')
        res = os.path.exists(os.path.join(output_dir, 'x0034_0B_x0020_0B.json'))
        self.assertTrue(res, True)
        cleanup_files(output_dir)

    def test_substructure_check_failing(self):
        synthons = ['c1ccccc1']
        fragmentA = get_mol('nsp13', 'x0176_0B', True)
        fragmentB = get_mol('nsp13', 'x0276_0B', True)
        filtered_synthons = merger.substructure_check(synthons, fragmentA, fragmentB)
        self.assertEqual(filtered_synthons, [])

    def test_substructure_check_passing(self):
        synthons = ['c1ccccc1']
        fragmentA = get_mol('nsp13', 'x0176_0B', True)
        fragmentB = get_mol('nsp13', 'x0438_0B', True)
        filtered_synthons = merger.substructure_check(synthons, fragmentA, fragmentB)
        self.assertTrue('c1ccccc1' in filtered_synthons)


if __name__ == "__main__":
    unittest.main()
