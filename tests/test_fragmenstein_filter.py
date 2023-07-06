"""Tests the Fragmenstein filter script"""

import os
import shutil
import unittest
from unittest.mock import patch

from filter.fragmenstein_filter import FragmensteinFilter, main
from utils.utils import get_mol, get_protein


def remove_files(dir):
    files = os.listdir(dir)
    files = [os.path.join(dir, f) for f in files]
    for file in files:
        try:
            shutil.rmtree(file)
        except NotADirectoryError:
            os.remove(file)


class TestFragmensteinFilter(unittest.TestCase):
    """Tests the fragmenstein filter functions"""

    def test_place_smiles(self):
        """Checks that molecules correctly pass and fail the filter"""
        frag_dir = os.path.join("tests", "test_Fragalysis")
        fragmentA_path = get_mol("Mpro", "x0107_0A", False, frag_dir)
        fragmentB_path = get_mol("Mpro", "x0678_0A", False, frag_dir)
        proteinA_path = get_protein("Mpro", "x0107_0A", False, frag_dir)
        proteinB_path = get_protein("Mpro", "x0678_0A", False, frag_dir)
        merge = "x107-0A-x0678-0A"
        name = "x107-0A-x0678-0A-123"
        smi = "NC(=O)CN1CCC2(C1)CC1(C2)OCCO1"
        synthon = "NC(=O)C[Xe]"
        working_dir = "tests/test_working/"
        output_dir = "tests/test_output/"

        filter = FragmensteinFilter(
            [smi],
            [synthon],
            fragmentA_path,
            fragmentB_path,
            proteinA_path,
            proteinB_path,
            merge,
            names=[name],
            work_pair_dir=working_dir,
            out_pair_dir=output_dir,
        )
        res = filter.filter_all()
        self.assertEqual(res[0][0], False)
        remove_files(working_dir)
        remove_files(output_dir)

    def test_main(self):
        frag_dir = os.path.join("tests", "test_Fragalysis")
        fragmentA_path = get_mol("Mpro", "x0107_0A", False, frag_dir)
        fragmentB_path = get_mol("Mpro", "x0678_0A", False, frag_dir)
        proteinA_path = get_protein("Mpro", "x0107_0A", False, frag_dir)
        proteinB_path = get_protein("Mpro", "x0678_0A", False, frag_dir)
        working_dir = "tests/test_working/"
        output_dir = "tests/test_output/"
        test_sdf = os.path.join("tests", "test_data", "fragmenstein_filter_mols.sdf")
        output_sdf = os.path.join("tests", "test_data", "test_fragmenstein_output.sdf")
        with patch(
            "sys.argv",
            [
                "filter/fragmenstein_filter.py",
                "-i",
                test_sdf,
                "-o",
                output_sdf,
                "-a",
                fragmentA_path,
                "-b",
                fragmentB_path,
                "-A",
                proteinA_path,
                "-B",
                proteinB_path,
                "-W",
                working_dir,
                "-O",
                output_dir,
                "-p",
                "x107-0A-x0678-0A"

            ]
        ):
            main()
        self.assertTrue(os.path.exists(output_sdf))
        os.remove(output_sdf)
        remove_files(working_dir)
        remove_files(output_dir)


if __name__ == "__main__":
    unittest.main()
