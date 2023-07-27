"""Perform R group expansion of a given merge"""

import argparse
import json

from fragment_network_merges.merge.find_merges import getFragmentNetworkSearcher


def perform_R_group_expansion(target, fragmentA, fragmentB, merges):
    fragmentNetworkSearcher = getFragmentNetworkSearcher()

    R_group_expansions = {}
    for merge in merges:
        expansions = fragmentNetworkSearcher.R_group_expansion(target,
                                                               merge,
                                                               fragmentA,
                                                               fragmentB)
        R_group_expansions[merge] = expansions
        print(len(expansions), 'expansions found for merge:', merge)
    return R_group_expansions


def perform_add_substituent(merge, substituent):
    fragmentNetworkSearcher = getFragmentNetworkSearcher()

    expansions = fragmentNetworkSearcher.add_known_substituent(merge,
                                                         substituent)
    return expansions


def main():
    parser = argparse.ArgumentParser(
        epilog="""
        Example
        python -m merge.R_group_expansions -t nsp13 -A x0034_0B -B x0176_0B -m c1ccccc1 c1ccc1
        OR
        python -m merge.R_group_expansions -M c1cccc1 -s [Xe]Cl
        """
    )
    parser.add_argument(
        "-t", "--target", help="the protein target (e.g. nsp13)", required=False
    )
    parser.add_argument(
        "-A",
        "--fragmentA_name",
        help="Fragalysis name for fragment A, e.g. x0034_0B",
        required=False
    )
    parser.add_argument(
        "-B",
        "--fragmentB_name",
        help="Fragalysis name for fragment B, e.g. x0034_0N",
        required=False
    )
    parser.add_argument(
        "-m",
        "--merges",
        nargs="+",
        help="list of merge SMILES for expansions",
        required=False
    )
    parser.add_argument(
        "-M",
        "--merge",
        help="single merge SMILES for expansion",
        required=False
    )
    parser.add_argument(
        "-s",
        "--substituent",
        help="SMILES of substituent to attempt expansion with (requires knowing attachment point (xenon) position) e.g. [Xe]Cl",
        required=False
    )
    parser.add_argument(
        "-o",
        "--output_fname",
        help="file to save expansions to",
        required=True
    )
    args = parser.parse_args()
    if args.substituent and args.merge:
        print('Expanding single SMILES')
        expansions = perform_add_substituent(args.merge, args.substituent)

    elif args.target and args.merges and args.fragmentA_name and args.fragmentB_name:
        print('Generating expansions for several SMILES')
        expansions = perform_R_group_expansion(args.target, args.fragmentA_name, args.fragmentB_name, args.merges)

    with open(args.output_fname, 'w') as f:
        json.dump(expansions, f)

    print('Expansions saved to file', args.output_fname)


if __name__ == "__main__":
    main()
