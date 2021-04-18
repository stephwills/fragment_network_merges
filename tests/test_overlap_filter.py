"""Tests the overlap filter script"""
import unittest
import sys
from rdkit import Chem
from rdkit.Chem import rdmolfiles
import numpy as np

# import scripts
sys.path.insert(1, '/home/sabsr3/xchem/fragment_network_merges/scripts')
from overlap_filter import OverlapFilter

# some test cases
suppl = Chem.SDMolSupplier('tests/test_embedded_mols.sdf')
mols = [x for x in suppl]
passing_mol = mols[0]
failing_mol = mols[1]
proteinA = rdmolfiles.MolFromPDBFile('tests/Mpro-x0107_0A_apo-desolv.pdb')
proteinB = rdmolfiles.MolFromPDBFile('tests/Mpro-x0678_0A_apo-desolv.pdb')

passing_case = OverlapFilter(passing_mol, proteinA, proteinB)
failing_case = OverlapFilter(failing_mol, proteinA, proteinB)

class TestOverlapFilter(unittest.TestCase):
    """Tests the OverlapFilter class"""

    def test_calc_distances(self):
        """Tests the function calculate distances correctly"""
        self.assertEqual((0.934984520123839, 0.9363730036063884), passing_case.calc_distances())
    
    def test_geometric_mean(self):
        """Tests the function calculates geometric mean"""
        mean = np.sqrt(0.934984520123839 * 0.9363730036063884)
        dists = passing_case.calc_distances()
        self.assertEqual(mean, passing_case.geometric_mean(dists[0], dists[1]))
    
    def test_filter(self):
        """Checks that molecules correctly pass and fail the filter"""
        self.assertEqual('pass', passing_case.filter())
        self.assertEqual('fail', failing_case.filter())

if __name__ == '__main__':
    unittest.main()
