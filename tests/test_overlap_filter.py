"""Tests the overlap filter script"""
import unittest
import numpy as np
from rdkit import Chem
from rdkit.Chem import rdmolfiles
from scripts.overlap_filter import calc_distances, geometric_mean, overlap_filter

# some test cases
suppl = Chem.SDMolSupplier('tests/test_embedded_mols.sdf')
mols = [x for x in suppl]
passing_mol = mols[0]
failing_mol = mols[1]
proteinA = rdmolfiles.MolFromPDBFile('tests/Mpro-x0107_0A_apo-desolv.pdb')
proteinB = rdmolfiles.MolFromPDBFile('tests/Mpro-x0678_0A_apo-desolv.pdb')

class TestOverlapFilter(unittest.TestCase):
    """Tests the overlap filter functions"""

    def test_calc_distances(self):
        """Tests the function calculate distances correctly"""
        self.assertEqual((0.934984520123839, 0.9363730036063884), calc_distances(passing_mol, proteinA, proteinB))
    
    def test_geometric_mean(self):
        """Tests the function calculates geometric mean"""
        mean = np.sqrt(0.934984520123839 * 0.9363730036063884)
        dists = calc_distances(passing_mol, proteinA, proteinB)
        self.assertEqual(mean, geometric_mean(dists[0], dists[1]))
    
    def test_filter(self):
        """Checks that molecules correctly pass and fail the filter"""
        self.assertEqual('pass', overlap_filter(passing_mol, proteinA, proteinB))
        self.assertEqual('fail', overlap_filter(failing_mol, proteinA, proteinB))

if __name__ == '__main__':
    unittest.main()
