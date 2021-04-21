"""Runs the preprocessing of fragments and merge generation"""

from joblib import Parallel, delayed
from scripts.find_merges import *
from scripts.preprocessing import *

# target = 'XXX'
# chain = 'XXX'
# fragments = ['XXX', 'XXX', 'XXX']

def preprocess_fragments(target, chain, fragment_names):
    """
    Preprocess the list of fragments and get all the possible fragment pairs
    to generate merges from.
    """
    # get all the smiles of the fragments
    fragment_smiles = [get_smiles(target, f, chain) for f in fragment_names]
    # check that the smiles exist as nodes in the network
    fragment_smiles, fragment_names = filter_for_nodes(fragment_smiles, fragment_names)
    # get all possible combinations of smiles
    smiles_pairs, name_pairs = get_combinations(fragment_smiles, fragment_names)
    smiles_pairs, name_pairs = check_fragment_pairs(smiles_pairs, name_pairs, target, chain)
    return smiles_pairs, name_pairs

# smiles_pairs, name_pairs = proprocess_fragments(target, chain, fragments)
# Parallel(n_jobs = 4)(delayed(get_expansions)(smiles_pair, name_pair) for smiles_pair, name_pair in zip(smiles_pairs, name_pairs))

fragment_names = ['x0104', 'x0107']
print(preprocess_fragments('Mpro', '0A', fragment_names))