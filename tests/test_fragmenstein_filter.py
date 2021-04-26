"""Tests the Fragmenstein filter script"""
import unittest
from rdkit import Chem
from scripts.fragmenstein_filter import *

passing_case = 'tests/test_dictionary.json'
failing_case = 'tests/test_dictionary2.json'

class TestFragmensteinFilter(unittest.TestCase):
    """Tests the fragmenstein filter functions"""

    def test_filter(self):
        """Checks that molecules correctly pass and fail the filter"""
        self.assertEqual(fragmenstein_filter(passing_case), 'pass')
        self.assertEqual(fragmenstein_filter(failing_case), 'fail')

if __name__ == '__main__':
    unittest.main()
