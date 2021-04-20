"""Runs the merge generation and filtering"""

import os
from joblib import Parallel, delayed
from rdkit import Chem
from rdkit.Chem import rdmolfiles
from scripts.find_merges import *
from scripts.descriptor_filter import *
from scripts.embedding_filter import *
from scripts.overlap_filter import *
from scripts.fragmenstein_filter import *
from scripts.interaction_fp_filter import *

def process_one_smi(name, smiles, fragmentA, fragmentB, synthon, proteinA, proteinB, output_directory):
    """
    Run all the filters on one smiles. If the smiles fails a filter, the function returns None.
    If it passes all the filters, it returns the smiles.

    :param name: the name of the merge (for naming output files - needs to be unique for each merge)
    :type name: string
    :param smiles: the smiles string of the merge
    :type smiles: string
    :param fragmentA: the filepath for fragment A mol file
    :type fragmentA: mol file
    :param fragmentB: the filepath for fragment B mol file
    :type fragmentB: mol file
    :param synthon: the smiles of the synthon used for the expansion
    :type synthon: string
    :proteinA: the filepath to the protein structure associated with fragment A
    :type proteinA: pdb file
    :proteinB: the filepath to the protein structure associated with fragment B
    :type proteinB: pdb file
    :param output_directory: directory for the fragmenstein output files
    :type output_directory: filepath

    :return: smiles (if all filters successful) or None
    :rtype: string or None
    """
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

results = Parallel(n_jobs = 4)(delayed(process_one_smi)(smi) for smi in list_of_smiles))
