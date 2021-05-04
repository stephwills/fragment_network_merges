"""Runs the filtering of the merges for one pair"""

import os
import argparse
from joblib import Parallel, delayed
from rdkit import Chem
from rdkit import RDLogger
from rdkit.Chem import rdmolfiles
from scripts.preprocessing import *
from scripts.descriptor_filter import *
from scripts.embedding_filter import *
from scripts.overlap_filter import *
from scripts.fragmenstein_filter import *
from scripts.interaction_fp_filter import *

RDLogger.DisableLog('rdApp.*')  # disable rdkit warning that occurs when mols can't be embedded

# get the file containing all the fragment pairs from the command line 
parser = argparse.ArgumentParser()
parser.add_argument('-f', '--merge_file', help='the json file containing the merges')
parser.add_argument('-m', '--merge', help='the name of the merge, e.g. x0107_0A_x0434_0A')
parser.add_argument('-a', '--fragment_A', help='fragment A mol file')
parser.add_argument('-b', '--fragment_B', help='fragment B mol file')
parser.add_argument('-p', '--protein_A', help='protein pdb file associated with fragment A')
parser.add_argument('-q', '--protein_B', help='protein pdb file associated with fragment B')
parser.add_argument('-o', '--output_directory', help='the directory to write the filtered files to')

# get the arguments
args = parser.parse_args()

# get the dictionary of merges - keys represent synthons, smiles 
merges_dict = open_json(args.merge_file)
synthons, smiles = get_merges(merges_dict)
num = range(len(smiles))  # used to give all the smiles an identifier

# get all the files needed for the function
fragmentA = args.fragment_A
fragmentB = args.fragment_B
proteinA = args.protein_A
proteinB = args.protein_B
output_directory = args.output_directory
merge_name = args.merge

def process_one_smi(num, smiles, synthon):
    """
    Run all the filters on one smiles. If the smiles fails a filter, the function returns None.
    If it passes all the filters, it returns the smiles.

    :param num: smiles are numbered to create an identifier (for naming output files etc.)
    :type num: int
    :param smiles: the smiles string of the merge
    :type smiles: string
    :param synthon: the smiles of the synthon used for the expansion
    :type synthon: string

    :return: smiles (if all filters successful) or None
    :rtype: string or None
    """
    # create unique 'name' for the molecule including the fragment merge and smiles number
    name = merge_name + '_' + str(num)
    # get molecules (some functions take in molecule rather than file)
    merge_mol = Chem.MolFromSmiles(smiles)
    fragmentA_mol = rdmolfiles.MolFromMolFile(fragmentA)
    fragmentB_mol = rdmolfiles.MolFromMolFile(fragmentB)
    synthon_mol = Chem.MolFromSmiles(synthon)
    proteinA_mol = rdmolfiles.MolFromPDBFile(proteinA)
    proteinB_mol = rdmolfiles.MolFromPDBFile(proteinB)
    # run the descriptor filter
    result = descriptor_filter(smiles)
    if result == 'fail':
        return None
    else:
        # run the embedding filter
        embedded, result = embedding_filter(merge_mol, fragmentA_mol, fragmentB_mol, synthon_mol)
        if result == 'fail':
            return None
        else:
            # run the overlap filter
            result = overlap_filter(embedded, proteinA_mol, proteinB_mol)
            if result == 'fail':
                return None
            else:
                # place with fragmenstein and run filter
                try:
                    json_fpath, placed_fpath = place_smiles(name, smiles, fragmentA, fragmentB, proteinA, output_directory)  # these are filenames
                    result = fragmenstein_filter(json_fpath)
                    if result == 'fail':
                        return None
                    else:
                        return smiles
                except:
                    return None
                # else:
                #     # run the interaction fp filter
                #     result = similarity_filter(placed_fpath, fragmentA, fragmentB, proteinA)  # these are filenames
                #     if result == 'fail':
                #         return None
                #     else:
                #         return smiles

# to test
# results = []
# for n, smi, syn in zip(num, smiles, synthons):
#     res = process_one_smi(n, smi, syn)
#     results.append(res)

results = Parallel(n_jobs = 8)(delayed(process_one_smi)(n, smi, syn) for n, smi, syn in zip(num, smiles, synthons))
filtered_results = [i for i in results if i]  # remove None values from list

filename = f'{output_directory}/{merge_name}_filtered.json'
with open(filename, 'w') as f:
    json.dump(filtered_results, f)
