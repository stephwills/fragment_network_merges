"""Tests the overlap filter script"""

import os
import unittest
import numpy as np

from rdkit import Chem

from filter.overlap_filter import OverlapFilter
from merge.preprocessing import get_protein


frag_dir = os.path.join('tests', 'test_Fragalysis')
suppl = Chem.SDMolSupplier(os.path.join('tests', 'test_data', 'overlap_filter_mols.sdf'))
test_mols = [x for x in suppl]
passing_mol, failing_mol = test_mols[0], test_mols[1]
proteinA = get_protein('Mpro', 'x0107_0A', frag_dir, True)
proteinB = get_protein('Mpro', 'x0678_0A', frag_dir, True)
proteinA_path = get_protein('Mpro', 'x0107_0A', frag_dir)
proteinB_path = get_protein('Mpro', 'x0678_0A', frag_dir)
filter = OverlapFilter(None, None, None, None, proteinA_path, proteinB_path, None, test_mols)

class TestOverlapFilter(unittest.TestCase):
    """Tests the overlap filter functions"""

    def test_calc_distances(self):
        """Tests the function calculate distances correctly"""
        distances = filter.calc_distances(passing_mol, proteinA, proteinB)
        actual_distances = (0.9268544843628844, 0.9311237700673226)
        self.assertEqual(distances, actual_distances)
    
    def test_geometric_mean(self):
        """Tests the function calculates geometric mean"""
        mean = np.sqrt(0.9268544843628844 * 0.9311237700673226)
        dists = filter.calc_distances(passing_mol, proteinA, proteinB)
        self.assertEqual(mean, filter.geometric_mean(dists[0], dists[1]))
    
    def test_filter_passing(self):
        """Checks that molecules correctly pass the filter"""
        passing_result = filter.filter_smi(passing_mol, proteinA, proteinB, 0.15)
        self.assertEqual(True, passing_result)

    def test_filter_failing(self):
        """Checks that molecules correctly fail the filter"""
        failing_result = filter.filter_smi(failing_mol, proteinA, proteinB, 0.15)
        self.assertEqual(False, failing_result)


if __name__ == '__main__':
    unittest.main()
