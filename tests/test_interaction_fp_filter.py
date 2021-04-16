import unittest
import sys
from rdkit import Chem

# import scripts
sys.path.insert(1, '/home/sabsr3/xchem/fragment_network_merges/scripts')
from interaction_fp_filter import InteractionFPFilter

class TestInteractionFPFilter(unittest.TestCase):

    def test_make_fp