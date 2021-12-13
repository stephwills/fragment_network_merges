"""Tests the Fragmenstein filter script"""

import os
import unittest
import shutil

from merge.preprocessing import get_mol, get_protein
from filter.fragmenstein_filter import FragmensteinFilter


passing_case = 'tests/test_dictionary.json'
failing_case = 'tests/test_dictionary2.json'
frag_dir = os.path.join('tests', 'test_Fragalysis')
fragmentA_path = get_mol('Mpro', 'x0107_0A', frag_dir)
fragmentB_path = get_mol('Mpro', 'x0678_0A', frag_dir)
proteinA_path = get_protein('Mpro', 'x0107_0A', frag_dir)
proteinB_path = get_protein('Mpro', 'x0678_0A', frag_dir)
merge = 'x107-0A-x0678-0A'
smi = 'NC(=O)CN1CCC2(C1)CC1(C2)OCCO1'
synthon = 'NC(=O)C[Xe]'
filter = FragmensteinFilter([smi], [synthon], fragmentA_path, fragmentB_path, proteinA_path, proteinB_path, merge)
working_dir = 'tests/test_working'
output_dir = 'tests/test_output'


def remove_files():
    if os.path.exists(os.path.join(working_dir, 'tempfiles')):
        shutil.rmtree(os.path.join(working_dir, 'tempfiles'))
    if os.path.exists(os.path.join(output_dir, 'fragmenstein')):
        shutil.rmtree(os.path.join(output_dir, 'fragmenstein'))


class TestFragmensteinFilter(unittest.TestCase):
    """Tests the fragmenstein filter functions"""

    def test_place_smiles(self):
        """Checks that molecules correctly pass and fail the filter"""
        filter._create_directories(working_dir, output_dir)
        filter.place_smiles(f'{merge}-4', smi, working_dir, output_dir, '2B')
        output_file = os.path.join(output_dir, 'fragmenstein', f'{merge}-4.minimised.mol')
        file_written = os.path.exists(output_file)
        self.assertEqual(file_written, True)
        remove_files()


if __name__ == '__main__':
    unittest.main()
