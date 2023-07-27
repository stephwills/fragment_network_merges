"""Tests the expansion filter script"""

import os
import unittest
from unittest.mock import patch

from fragment_network_merges.filter import ExpansionFilter, main
from fragment_network_merges.merge import get_mol

fragalysis_dir = "tests/test_Fragalysis"
fragmentA = get_mol("nsp13", "x0276_0B", True)
fragmentB = get_mol("nsp13", "x0034_0B", True)
test_sdf = os.path.join("tests", "test_data", "embedding_filter_mols.sdf")
output_sdf = os.path.join("tests", "test_data", "test_embedding_output.sdf")
test_fragmentA = get_mol("nsp13", "x0034_0B", False, fragalysis_dir)
test_fragmentB = get_mol("nsp13", "x0176_0B", False, fragalysis_dir)


class TestExpansionFilter(unittest.TestCase):
    """Tests the expansion filter function"""

    def test_expansion_filter_passing(self):
        """Tests the filter passes molecules correctly"""
        smi = "Fc1ccc(CCNc2ccc(F)cc2)cc1"
        synthon = "Fc1ccccc1[Xe]"
        filter = ExpansionFilter()
        passing_case = filter.filter_smi(smi, synthon, fragmentA, fragmentB, 3)
        self.assertEqual(passing_case, True)

    def test_expansion_filter_valence(self):
        """Tests the filter for molecules that have create valence errors"""
        smis = [
            "CCC(Cc1ccco1)NCCn1cnc(C)c1C",
            "Cn1c(-c2ccco2)nn(CN2CCc3ccccc3C2)c1=S",
            "Cc1ncn(C(C)C(=O)NCc2ccco2)c1C",
            "Cc1[nH]c(-c2ccco2)nc1CN1CCNCC1",
            "Cc1ncn(CC(=O)N(C)Cc2ccco2)c1C",
            "Cc1ncn(CCNC(=O)N(Cc2ccncc2)CC(C)C)c1C",
            "Cc1ncn(CCNC(=O)Cc2ccncc2)c1C",
            "CCC(NC(=O)Cn1cnc(C)c1C)c1ccncc1",
        ]
        syns = [
            "[Xe]c1ccco1",
            "[Xe]c1ccco1",
            "[Xe]c1ccco1",
            "[Xe]c1ccco1",
            "[Xe]c1ccco1",
            "[Xe]c1ccncc1",
            "[Xe]c1ccncc1",
            "[Xe]c1ccncc1",
        ]
        fragmentA = get_mol("PGN_RS02895PGA", "x0051_0A", True)
        fragmentB = get_mol("PGN_RS02895PGA", "x0087_0A", True)
        filter = ExpansionFilter()
        results = [
            filter.filter_smi(smi, syn, fragmentA, fragmentB, 3)
            for smi, syn in zip(smis, syns)
        ]
        test_results = [True] * 8
        self.assertEqual(results, test_results)

    def test_main(self):
        """Check that main executes correctly and produces output file"""
        with patch(
            "sys.argv",
            [
                "filter/expansion_filter.py",
                "-i",
                test_sdf,
                "-o",
                output_sdf,
                "-A",
                test_fragmentA,
                "-B",
                test_fragmentB,
                "--min_atoms",
                "2",
            ],
        ):
            main()
        self.assertTrue(os.path.exists(output_sdf))
        os.remove(output_sdf)


if __name__ == "__main__":
    unittest.main()
