"""Runs the preprocessing of fragments and merge generation"""

import argparse

from scripts.preprocessing import get_smiles, check_fragment_pairs
from scripts.find_merges import filter_for_nodes, get_combinations, get_expansions


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

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-t', '--target', help='the protein target (e.g. nsp13)')
    parser.add_argument('-f', '--fragments', nargs='+', help='the list of fragments to merge')
    parser.add_argument('-o', '--output_directory', help='the directory to write the merge files to')
    args = parser.parse_args()

    # get all fragment pairs and check they exist in the network
    smiles_pairs, name_pairs = preprocess_fragments(args.target, args.fragments)

    # run the database query; this will create json files containing the merges for each pair
    # in the output directory
    for smiles_pair, name_pair in zip(smiles_pairs, name_pairs):
        get_expansions(smiles_pair, name_pair, args.target, args.output_directory)

if __name__ == "__main__":
    main()
