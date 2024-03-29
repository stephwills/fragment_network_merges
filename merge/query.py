"""Runs the preprocessing of fragments and merge generation"""

import argparse

from merge.config_merge import config_merge
from merge.find_merges import getFragmentNetworkSearcher
from merge.preprocessing import check_fragment_pairs, check_merges_run, check_too_similar
from utils.utils import get_smiles


def preprocess_fragments(target, fragment_names, output_dir=None, remove_similar=False):
    """
    Preprocess the list of fragments and get all the possible fragment pairs
    to generate merges from.
    """
    fragmentNetworkSearcher = getFragmentNetworkSearcher()

    # get all the smiles of the fragments
    fragment_smiles = [get_smiles(target, f) for f in fragment_names]

    # check that the smiles exist as nodes in the network
    fragment_smiles, fragment_names = fragmentNetworkSearcher.filter_for_nodes(
        fragment_smiles, fragment_names
    )

    # get all possible combinations of smiles
    smiles_pairs, name_pairs = fragmentNetworkSearcher.get_combinations(
        fragment_smiles, fragment_names
    )
    smiles_pairs, name_pairs = check_fragment_pairs(
        smiles_pairs, name_pairs, target
    )  # check distance between pair
    if remove_similar:
        smiles_pairs, name_pairs = check_too_similar(
            smiles_pairs,
            name_pairs,
            target
        )
    smiles_pairs, name_pairs = check_merges_run(
        smiles_pairs, name_pairs, output_dir
    )  # check if merges run already

    return smiles_pairs, name_pairs


def main():
    parser = argparse.ArgumentParser(
        epilog="""
    Example
    python -m merge.query -t nsp13 -f x0034_0B x0176_0B x0212_0B -o data/example_folder"""
    )
    parser.add_argument(
        "-t", "--target", help="the protein target (e.g. nsp13)", required=True
    )
    parser.add_argument(
        "-f",
        "--fragments",
        nargs="+",
        help="the list of fragments to merge. E.g. x0032_0A x0034_0B",
        required=True,
    )
    parser.add_argument(
        "-o",
        "--output_dir",
        help="the directory to write the merge files to",
        required=True,
    )
    parser.add_argument(
        "-w",
        "--working_dir",
        help="the directory where intermediate results will be computed",
        required=False,
    )
    parser.add_argument(
        "-s",
        "--remove_similar_fragments",
        help="remove fragments that are too similar from merging query",
        default=False
    )

    args = parser.parse_args()

    # check if working dir set
    if args.working_dir:
        config_merge.WORKING_DIR = args.working_dir

    # get all fragment pairs and check they exist in the network
    smiles_pairs, name_pairs = preprocess_fragments(
        args.target, args.fragments, args.output_dir, args.remove_similar_fragments
    )
    print("Fragments have been processed!")

    # run the database query; this will create json files containing the merges for each pair
    # in the output directory
    fragmentNetworkSearcher = getFragmentNetworkSearcher()
    fragmentNetworkSearcher.expand_fragmentA(
        name_pairs, args.target, output_dir=args.output_dir
    )


if __name__ == "__main__":
    main()
