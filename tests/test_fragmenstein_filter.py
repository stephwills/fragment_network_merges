import unittest
import sys
from rdkit import Chem

# import scripts
sys.path.insert(1, '/home/sabsr3/xchem/fragment_network_merges/scripts')
from fragmenstein_filter import FragmensteinFilter

passing_case = FragmensteinFilter('tests/test_dictionary.json')
failing_case = FragmensteinFilter('tests/test_dictionary2.json')

class TestFragmensteinFilter(unittest.TestCase):

    def test_filter(self):
        self.assertEqual(passing_case.filter(), 'pass')
        self.assertEqual(failing_case.filter(), 'fail')

if __name__ == '__main__':
    unittest.main()
