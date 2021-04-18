"""Tests the Fragmenstein filter script"""
import unittest
from rdkit import Chem
from scripts.fragmenstein_filter import FragmensteinFilter

passing_case = FragmensteinFilter('tests/test_dictionary.json')
failing_case = FragmensteinFilter('tests/test_dictionary2.json')

class TestFragmensteinFilter(unittest.TestCase):
    """Tests the FragmensteinFilter class"""

    def test_filter(self):
        """Checks that molecules correctly pass and fail the filter"""
        self.assertEqual(passing_case.filter(), 'pass')
        self.assertEqual(failing_case.filter(), 'fail')

if __name__ == '__main__':
    unittest.main()
