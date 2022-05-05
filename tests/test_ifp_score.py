"""Tests the interaction fingeprint scoring script"""

import os
import unittest

from filter.ifp_score import IfpScore
from utils.utils import get_mol


class TestIfpScore(unittest.TestCase):
    """Tests the IFP scoring function"""

    def test_scoring(self):
        """
        Tests the function scores molecules correctly
        """
        frag_dir = os.path.join("tests", "test_Fragalysis")
        fragmentA = get_mol("nsp13", "x0034_0B", False, frag_dir)
        fragmentB = get_mol("nsp13", "x0212_0B", False, frag_dir)
        names = [
            "x0034_0B_x0212_0B_16954",
            "x0034_0B_x0212_0B_17804",
            "x0034_0B_x0212_0B_18254",
        ]
        mol_dir = os.path.join("tests", "test_data", "for_scoring")
        mol_fnames = [n.replace("_", "-") for n in names]
        mol_files = [os.path.join(mol_dir, f"{f}.minimised.mol") for f in mol_fnames]
        apo_files = [os.path.join(mol_dir, f"{f}.holo_minimised_nolig.pdb") for f in mol_fnames]
        scorer = IfpScore(
            fragmentA=fragmentA,
            fragmentB=fragmentB,
            mol_files=mol_files,
            apo_files=apo_files,
        )

        scores = scorer.score_all()
        actual_scores = [0.38, 0.62, 0.58]
        rounded_scores = [round(s, 2) for s in scores]
        self.assertEqual(actual_scores, rounded_scores)


if __name__ == "__main__":
    unittest.main()
