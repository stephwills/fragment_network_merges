"""Tests the interaction fingeprint scoring script"""

import os
import unittest

from filter.plip_ifp_score import PlipIfpScore
from utils.utils import get_mol


class TestPlipIfpScore(unittest.TestCase):
    """Tests the PLIP IFP scoring function"""

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
        holo_files = [os.path.join(mol_dir, f"{f}.holo_minimised.pdb") for f in mol_fnames]
        work_dir = mol_dir
        output_dir = os.path.join("tests", "test_output")
        scorer = PlipIfpScore(
            names=names,
            fragmentA=fragmentA,
            fragmentB=fragmentB,
            mol_files=mol_files,
            apo_files=apo_files,
            holo_files=holo_files,
            work_pair_dir=work_dir,
            out_pair_dir=output_dir
        )

        scores = scorer.score_all()
        rounded_scores = [0.5, 0.6, 0.6]
        self.assertEqual(rounded_scores, [round(s, 1) for s in scores])


if __name__ == "__main__":
    unittest.main()
