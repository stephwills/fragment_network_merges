"""Utils functions for similarity search"""

import json
import os

from utils.utils import get_smiles, load_json


def process_sim_search_data(
    target: str,
    raw_file: str,
    fragments: list,
    output_dir: str,
    tvThreshold: float = 0.4,
    returnData=False,
):
    """
    Process similarity search data so that it can be input into the filtering pipeline. Creates separate json files
    for each pair of fragments containing a list of the SMILES of the merges. If specified, a Tversky threshold
    is applied as a filter.
    """
    # similarity search uses the SMILES of the fragments - create dict to retrieve fragment names for filtering
    fragment_dict = {}

    for fragment in fragments:
        smi = get_smiles(target, fragment)
        fragment_dict[smi] = fragment

    tverskys = []
    ids = []
    smiles = []
    fAs = []
    fBs = []
    pairs = []

    raw_data = load_json(raw_file)
    for pair in raw_data:
        frags = pair.split(",")
        frags.sort()
        pair_smiles = []
        smiA, smiB = frags[0], frags[1]
        fA, fB = fragment_dict[smiA], fragment_dict[smiB]
        output_fname = os.path.join(output_dir, target, f"{fA}_{fB}.json")
        for lst in raw_data[pair]:
            tv = lst[0]
            if tv >= tvThreshold:
                pair_smiles.append(lst[2])
                # record other data
                tverskys.append(tv)
                ids.append(lst[1])
                smiles.append(lst[2])
                pair_smiles.append(lst[2])
                fAs.append(fA)
                fBs.append(fB)
                pairs.append(f"{fA}_{fB}")

        with open(output_fname, "w") as f:
            json.dump(list(set(pair_smiles)), f)

    print(f"{len(smiles)} smiles found for target {target}")
    print(f"{len(set(smiles))} unique smiles found for target {target}")

    if returnData:
        return pairs, fAs, fBs, smiles, tverskys, ids

target = 'SDCBPA'
raw_file = '/home/swills/Oxford/data/fragment_merging/SDCBPA/similarity_search/raw_input/database.json'

fragments = ['x0230', 'x0441', 'x013_4', 'x013_2', 'x0456', 'x0251', 'x013_A', 'x013_3', 'x0144', 'x013_1', 'x0430', 'x013_5']
output_dir = '/home/swills/Oxford/data/fragment_merging/SDCBPA/similarity_search/test'
print(process_sim_search_data(target,raw_file, fragments,output_dir))