"""Tests the expansion filter script"""

import unittest

from filter.expansion_filter import ExpansionFilter
from merge.preprocessing import get_mol

fragmentA_path = get_mol("nsp13", "x0276_0B", False)
fragmentB_path = get_mol("nsp13", "x0034_0B", False)
fragmentA = get_mol("nsp13", "x0276_0B", True)
fragmentB = get_mol("nsp13", "x0034_0B", True)


class TestExpansionFilter(unittest.TestCase):
    """Tests the expansion filter function"""

    def test_expansion_filter_failing(self):
        """
        Tests the filter fails where most of the synthon is also in fragment A
        (so looks like an expansion).
        """
        smi = "CC(C)NCC(C)(C)CN(C)C(=O)c1ccccc1F"
        synthon = "Fc1ccccc1[Xe]"
        filter = ExpansionFilter([smi], [synthon], fragmentA_path, fragmentB_path)
        failing_case = filter.filter_smi(smi, synthon, fragmentA, fragmentB)
        self.assertEqual(failing_case, False)

    def test_expansion_filter_passing(self):
        """Tests the filter passes molecules correctly"""
        smi = "Fc1ccc(CCNc2ccc(F)cc2)cc1"
        synthon = "Fc1ccccc1[Xe]"
        filter = ExpansionFilter([smi], [synthon], fragmentA_path, fragmentB_path)
        failing_case = filter.filter_smi(smi, synthon, fragmentA, fragmentB)
        self.assertEqual(failing_case, True)

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
        fA_path = get_mol("PGN_RS02895PGA", "x0051_0A", False)
        fB_path = get_mol("PGN_RS02895PGA", "x0087_0A", False)
        filter = ExpansionFilter(smis, syns, fA_path, fB_path)
        results, _ = filter.filter_all(cpus=1)
        test_results = [True] * 8
        self.assertEqual(results, test_results)


if __name__ == "__main__":
    unittest.main()
