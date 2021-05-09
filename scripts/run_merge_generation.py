"""Runs the preprocessing of fragments and merge generation"""

from joblib import Parallel, delayed
from scripts.find_merges import *
from scripts.preprocessing import *

def preprocess_fragments(target, fragment_names):
    """
    Preprocess the list of fragments and get all the possible fragment pairs
    to generate merges from.
    """
    # get all the smiles of the fragments
    fragment_smiles = [get_smiles(target, f) for f in fragment_names]
    # check that the smiles exist as nodes in the network
    fragment_smiles, fragment_names = filter_for_nodes(fragment_smiles, fragment_names)
    # get all possible combinations of smiles
    smiles_pairs, name_pairs = get_combinations(fragment_smiles, fragment_names)
    smiles_pairs, name_pairs = check_fragment_pairs(smiles_pairs, name_pairs, target)
    return smiles_pairs, name_pairs

target = 'nsp13'
frags = ['x0034_0B', 'x0176_0B', 'x0183_0B', 'x0208_0A', 'x0212_0B', 'x0246_0B', 'x0276_0B', 'x0283_0B', 'x0311_0B', 'x0438_0B']

smiles_pairs, name_pairs = preprocess_fragments(target, frags)

for smiles_pair, name_pair in zip(smiles_pairs, name_pairs):
    get_expansions(smiles_pair, name_pair, target)
