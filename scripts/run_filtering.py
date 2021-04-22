"""Runs the filtering of the merges for one pair"""

import os
import argparse
from joblib import Parallel, delayed
from rdkit import Chem
from rdkit.Chem import rdmolfiles
from scripts.descriptor_filter import *
from scripts.embedding_filter import *
from scripts.overlap_filter import *
from scripts.fragmenstein_filter import *
from scripts.interaction_fp_filter import *

# get the file containing all the fragment pairs from the command line 
parser = argparse.ArgumentParser()
parser.add_argument('-m', '--merge_file', help='the json file containing the merges')
parser.add_argument('-a', '--fragment_A_file', help='the mol file for fragment A')
parser.add_argument('-b', '--fragment_B_file', help='the mol file for fragment B')
parser.add_argument('-p', '--protein_A_file', help='the pdb file containing the protein associated with fragment A')
parser.add_argument('-q', '--protein_B_file', help='the pdb file containing the protein associated with fragment B')
parser.add_argument('-o', '--output_directory', help='the directory to write the fragmenstein files to')

# get the arguments
args = parser.parse_args()

# get the dictionary of merges - keys represent synthons, smiles 
merges_dict = open_json(args.merge_file)
synthons, smiles = get_merges(merges_dict)
num = range(len(smiles))  # used to give all the smiles an identifier

# get all the files needed for the function
fragmentA = args.fragment_A_file
fragmentB = args.fragment_B_file
proteinA = args.protein_A_file
proteinB = args.protein_B_file
output_directory = args.output_directory

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
    merge_name = args.merge_file.replace('.json', '').replace('data/', '')
    name = merge_name + '_' + str(num)
    # get molecules (some functions take in molecule rather than file)
    merge_mol = Chem.MolFromSmiles(smiles)
    fragmentA_mol = rdmolfiles.MolFromMolFile(fragmentA)
    fragmentB_mol = rdmolfiles.MolFromMolFile(fragmentB)
    synthon_mol = Chem.MolFromSmiles(synthon)
    proteinA_mol = rdmolfiles.MolFromPDBFile(proteinA)
    proteinA_mol = rdmolfiles.MolFromPDBFile(proteinB)
    # run the descriptor filter
    result = descriptor_filter(smiles)
    if result == 'fail':
        return None
    # run the embedding filter
    embedded, result = embedding_filter(merge_mol, fragmentA_mol, fragmentB_mol, synthon_mol)
    if result == 'fail':
        return None
    # run the overlap filter
    result = overlap_filter(embedded, proteinA_mol, proteinB_mol)
    if result == 'fail':
        return None
    # place with fragmenstein and run filter
    place_smiles(name, smiles, fragmentA, fragmentB, proteinA, output_directory)  # these are filenames
    json_fname = name + '.minimised.json'
    json_fpath = os.path.join(output_directory, name, json_fname) 
    result = fragmenstein_filter(json_fname)
    if result == 'fail':
        return None
    # run the interaction fp filter
    placed_fname = name + '.minimised.mol'
    placed_path = os.path.join(output_directory, name, placed_fname)
    result = similarity_filter(placed_fname, fragmentA, fragmentB, proteinA)  # these are filenames
    if result == 'fail':
        return None
    else:
        return smiles

results = Parallel(n_jobs = 4)(delayed(process_one_smiles)(n, smi, syn) for n, smi, syn in zip(num, smiles, synthons))