"""Runs the preprocessing of fragments and merge generation"""

import argparse
import os

from scripts.config import config
from scripts.preprocessing import get_smiles, check_fragment_pairs
from scripts.find_merges import  getFragmentNetworkSearcher


def preprocess_fragments(target, fragment_names):
    """
    Preprocess the list of fragments and get all the possible fragment pairs
    to generate merges from.
    """
    fragmentNetworkSearcher = getFragmentNetworkSearcher()
    # get all the smiles of the fragments
    fragment_smiles = [get_smiles(target, f) for f in fragment_names]

    # check that the smiles exist as nodes in the network
    fragment_smiles, fragment_names = fragmentNetworkSearcher.filter_for_nodes(fragment_smiles, fragment_names)

    # get all possible combinations of smiles
    smiles_pairs, name_pairs = fragmentNetworkSearcher.get_combinations(fragment_smiles, fragment_names)
    smiles_pairs, name_pairs = check_fragment_pairs(smiles_pairs, name_pairs, target)
    return smiles_pairs, name_pairs

def main():
    parser = argparse.ArgumentParser(epilog='''
    Example
    python -m scripts.run_database_query -t nsp13 -f x0034_0B x0176_0B x0212_0B -o data/example_folder''')
    parser.add_argument('-t', '--target', help='the protein target (e.g. nsp13)', required=True)
    parser.add_argument('-f', '--fragments', nargs='+', help='the list of fragments to merge. E.g. x0032_0A x0034_0B', required=True)
    parser.add_argument('-o', '--output_directory', help='the directory to write the merge files to', required=True)
    parser.add_argument('-w', '--wdir', help='the directory where intermediate results will be computed', required=False)

    args = parser.parse_args()

    if os.path.isdir(args.target):
        config.FRAGALYSIS_DATA_DIR = os.path.split( args.target )[0]

    config.WORKING_DIR = args.wdir

    # get all fragment pairs and check they exist in the network
    # TODO: Ensure that already computed results are not recomputed.
    smiles_pairs, name_pairs = preprocess_fragments(args.target, args.fragments)
    print("Fragments have been processed!")

    fragmentNetworkSearcher = getFragmentNetworkSearcher()

    # run the database query; this will create json files containing the merges for each pair
    # in the output directory
    for smiles_pair, name_pair in zip(smiles_pairs, name_pairs): #TODO: We may perform several parallel queries, need to talk to Tim
        fragmentNetworkSearcher.get_expansions(smiles_pair, name_pair, args.target, args.output_directory)

if __name__ == "__main__":
    main()
