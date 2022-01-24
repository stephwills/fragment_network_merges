"""Tests the interaction fingeprint scoring script"""

import os
import unittest

from merge.preprocessing import get_mol, get_protein
from filter.sucos_score import SuCOSScore


smis = ['smi1', 'smi2', 'smi3']
synthons = ['synthon1', 'synthon2', 'synthon3']
frag_dir = os.path.join('tests', 'test_Fragalysis')
fragmentA = get_mol('nsp13', 'x0034_0B', frag_dir)
fragmentB = get_mol('nsp13', 'x0212_0B', frag_dir)
proteinA = get_protein('nsp13', 'x0034_0B', frag_dir)
proteinB = get_protein('nsp13', 'x0212_0B', frag_dir)
merge = 'x0034_0B_x0212_0B'
mols = ['mol1', 'mol2', 'mol3']
names = ['x0034_0B_x0212_0B_16954', 'x0034_0B_x0212_0B_17804', 'x0034_0B_x0212_0B_18254']
mol_dir = os.path.join('tests', 'test_data', 'for_scoring')
mol_fnames = [n.replace('_', '-') for n in names]
mol_files = [os.path.join(mol_dir, f"{f}.minimised.mol") for f in mol_fnames]
holo_files = [os.path.join(mol_dir, f"{f}.holo_minimised.pdb") for f in mol_fnames]
apo_files = [os.path.join(mol_dir, f"{f}.holo_minimised_nolig.pdb") for f in mol_fnames]
score = SuCOSScore(smis, synthons, fragmentA, fragmentB, proteinA, proteinB, merge, mols, names, mol_files,
                 holo_files, apo_files)


class TestSuCOSScore(unittest.TestCase):
    """Tests the SuCOS scoring function"""

    def test_scoring(self):
        """
        Tests the function scores molecules correctly
        """
        scores = score.score_all()
        actual_scores = [0.53, 0.71, 0.67]
        rounded_scores = [round(s, 2) for s in scores]
        self.assertEqual(actual_scores, rounded_scores)


if __name__ == '__main__':
    unittest.main()
