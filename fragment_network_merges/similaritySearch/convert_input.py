"""Convert input to be compatible with filtering pipeline"""
from argparse import ArgumentParser
import json
import os

from tqdm import tqdm
from fragment_network_merges.utils.utils import get_smiles, load_json


def process_input(fragments, target, file, output_dir, tv_threshold=0.4, return_data=False):
    """
    Process the output file from the similarity search so it can be used as input into the filtering pipeline.

    :param fragments: list of fragment names
    :type fragments: list
    :param target: name of the target, e.g. Mpro
    :type target: str
    :param file: name of the file that came out of the similarity search
    :type file: str
    :param output_dir: path to the dir to write the processed files for input into the filtering pipeline
    :type output_dir: str
    :param tv_threshold: mean Tversky threshold to select compounds
    :type tv_threshold: float
    """
    fragment_dict = {}
    for fragment in fragments:
        smi = get_smiles(target, fragment)
        fragment_dict[smi] = fragment

    tverskys = []
    ids = []
    smiles = []
    pairs = []

    data = load_json(file)

    for pair in tqdm(data):
        frags = pair.split(',')
        frags.sort()
        pair_smiles = []
        smiA, smiB = frags[0], frags[1]
        fA, fB = fragment_dict[smiA], fragment_dict[smiB]
        output_fname = os.path.join(output_dir, f"{fA}_{fB}.json")

        for lst in data[pair]:
            tv = lst[0]
            if tv >= tv_threshold:
                tverskys.append(tv)
                ids.append(lst[1])
                smiles.append(lst[2])
                pair_smiles.append(lst[2])
                pairs.append(f"{fA}_{fB}")

        with open(output_fname, 'w') as f:
            json.dump(list(set(pair_smiles)), f)

    print(f'{len(smiles)} smiles found for target {target}')
    print(f'{len(set(smiles))} unique smiles found for target {target}')

    if return_data:
        return smiles, ids, tverskys, pairs


def main():
    parser = ArgumentParser()
    parser.add_argument('-F', '--fragments', nargs='+', help='list of fragments for target')
    parser.add_argument('-t', '--target')
    parser.add_argument('-f', '--sim_search_file')
    parser.add_argument('-o', '--output_dir', help='where to write processed files')
    parser.add_argument('-T', '--tv_threshold', type=float, default=0.4, help='mean Tversky threshold for filtering compounds')
    args = parser.parse_args()

    process_input(args.fragments,
                  args.target,
                  args.sim_search_file,
                  args.output_dir,
                  args.tv_threshold)


if __name__ == "__main__":
    main()
