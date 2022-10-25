"""Tests the interaction fingeprint scoring script"""

import os
import unittest

from filter.elaboratability_score import ElaboratabilityScore
from utils.utils import get_mol
from rdkit import Chem


def delete_files(dir):
    for intermed_file in os.listdir(dir):
        if 'elab' in intermed_file:
            if os.path.exists(os.path.join(dir, intermed_file)):
                os.remove(os.path.join(dir, intermed_file))


class TestElaboratabilityScore(unittest.TestCase):
    """Tests the PLIP IFP scoring function"""

    def test_scoring(self):
        """
        Tests the function scores molecules correctly
        """
        names = [
            "x0034_0B_x0212_0B_16954",
            "x0034_0B_x0212_0B_17804",
            "x0034_0B_x0212_0B_18254",
        ]
        mol_dir = os.path.join("tests", "test_data", "for_scoring")
        mol_fnames = [n.replace("_", "-") for n in names]
        mol_files = [os.path.join(mol_dir, f"{f}.minimised.mol") for f in mol_fnames]
        mols = [Chem.MolFromMolFile(m) for m in mol_files]
        apo_files = [os.path.join(mol_dir, f"{f}.holo_minimised_nolig.pdb") for f in mol_fnames]
        work_dir = mol_dir
        output_dir = os.path.join("tests", "test_output")
        scorer = ElaboratabilityScore(
            names=names,
            mol_files=mol_files,
            apo_files=apo_files,
            mols=mols,
            work_pair_dir=work_dir,
            out_pair_dir=output_dir
        )

        scores = scorer.score_all()
        self.assertEqual(scores, [3,2,3])
        delete_files(mol_dir)


if __name__ == "__main__":
    unittest.main()
