"""Tests the interaction fingeprint scoring script"""

import os
import shutil
import unittest
from unittest.mock import patch

from filter.plip_ifp_score import PlipIfpScore, main
from utils.utils import get_mol


def delete_files(dir):
    for intermed_file in os.listdir(dir):
        if os.path.isdir(os.path.join(os.path.join(dir, intermed_file))):
            shutil.rmtree(os.path.join(os.path.join(dir, intermed_file)))
        for file_sub in ["interactions", "fA", "fB"]:
            if file_sub in intermed_file:
                if os.path.exists(os.path.join(dir, intermed_file)):
                    os.remove(os.path.join(dir, intermed_file))


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
        apo_files = [
            os.path.join(mol_dir, f"{f}.holo_minimised_nolig.pdb") for f in mol_fnames
        ]
        holo_files = [
            os.path.join(mol_dir, f"{f}.holo_minimised.pdb") for f in mol_fnames
        ]
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
            out_pair_dir=output_dir,
        )

        scores = scorer.score_all()
        rounded_scores = [0.5, 0.6, 0.6]
        self.assertEqual(rounded_scores, [round(s, 1) for s in scores])
        delete_files(mol_dir)

    def test_main(self):
        """Check that main executes correctly and produces output file"""
        input_sdf = os.path.join(
            "tests", "test_data", "for_scoring", "input_scoring.sdf"
        )
        output_sdf = os.path.join(
            "tests", "test_data", "for_scoring", "output_scoring.sdf"
        )
        frag_dir = os.path.join("tests", "test_Fragalysis")
        fragmentA = get_mol("nsp13", "x0034_0B", False, frag_dir)
        fragmentB = get_mol("nsp13", "x0212_0B", False, frag_dir)
        dir = os.path.join("tests", "test_data", "for_scoring")
        with patch(
            "sys.argv",
            [
                "filter/plip_ifp_score.py",
                "-i",
                input_sdf,
                "-o",
                output_sdf,
                "-a",
                fragmentA,
                "-b",
                fragmentB,
                "-t",
                "0.5",
                "-W",
                dir,
                "-O",
                dir,
            ],
        ):
            main()
        self.assertTrue(os.path.exists(output_sdf))
        os.remove(output_sdf)
        delete_files(os.path.join("tests", "test_data", "for_scoring"))


if __name__ == "__main__":
    unittest.main()
