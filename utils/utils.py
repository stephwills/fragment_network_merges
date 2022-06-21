"""
Utils functions for both pipelines.
"""

import json
import os

import numpy as np
from merge.config_merge import config_merge
from rdkit import Chem
from rdkit.Chem import rdmolfiles


def load_json(fname):
    """
    Function to open json file
    """
    with open(fname, "r") as f:
        data = json.load(f)
    return data


def get_smiles(target, fragment, fragalysis_dir=config_merge.FRAGALYSIS_DATA_DIR):
    """
    Function to get the SMILES for each fragment.
    File paths are like this: TARGET/aligned/TARGET-FRAGMENT_CHAIN
    e.g. Mpro/aligned/Mpro-x0464_0A

    :param target: the protein target, e.g. 'Mpro'
    :type target: string
    :param fragment: the fragment i.d., e.g. 'x0107_0A'
    :type fragment: string
    :param fragalysis_dir: the directory containing fragalysis data; within should be target folder, e.g. 'Mpro'
    :type fragalysis_dir: string

    :return: smiles string
    :rtype: string
    """
    dir = os.path.join(fragalysis_dir, target, "aligned")
    fname_part = f"{target}-{fragment}"  # the first part of the filenames/folder names
    path = os.path.join(dir, fname_part)
    smiles_path = os.path.join(path, f"{fname_part}_smiles.txt")

    try:
        with open(smiles_path) as smiles_file:  # open the file to get smiles
            smiles = smiles_file.read()
    except OSError as e:
        print("SMILES file cannot be found for that fragment.")
        print("Ensure files are saved correctly in Fragalysis format.")
        print(e)

    # smiles does not necessarily match what is in the network
    # get canonical smiles by converting to mol and converting back to smiles
    mol = Chem.MolFromSmiles(smiles)
    smiles = Chem.MolToSmiles(mol)
    return smiles


def get_mol(
    target, fragment, return_mol=False, fragalysis_dir=config_merge.FRAGALYSIS_DATA_DIR
):
    """
    Function to get the mol for each fragment.
    File paths are like this: TARGET/aligned/TARGET-FRAGMENT_CHAIN
    e.g. Mpro/aligned/Mpro-x0464_0A

    :param target: the protein target, e.g. 'Mpro'
    :type target: string
    :param fragment: the fragment i.d., e.g. 'x0107'
    :type fragment: string
    :param return_mol: whether to return RDKit molecule (True) or path (False)
    :type return_mol: bool
    :param fragalysis_dir: the directory containing fragalysis data; within should be target folder, e.g. 'Mpro'
    :type fragalysis_dir: string

    :return: smiles string
    """
    dir = os.path.join(fragalysis_dir, target, "aligned")
    fname_part = f"{target}-{fragment}"  # the first part of the filenames/folder names
    path = os.path.join(dir, fname_part)
    mol_path = os.path.join(path, f"{fname_part}.mol")

    if return_mol:
        try:
            mol = rdmolfiles.MolFromMolFile(mol_path)
            return mol
        except OSError as e:
            print("Mol file cannot be found for that fragment.")
            print("Ensure files are saved correctly in Fragalysis format.")
            print(e)
    else:
        return mol_path


def get_protein(
    target, fragment, return_mol=False, fragalysis_dir=config_merge.FRAGALYSIS_DATA_DIR
):
    """
    Function to get the mol for each fragment.
    File paths are like this: TARGET/aligned/TARGET-FRAGMENT_CHAIN
    e.g. Mpro/aligned/Mpro-x0464_0A

    :param target: the protein target, e.g. 'Mpro'
    :type target: string
    :param fragment: the fragment i.d., e.g. 'x0107'
    :type fragment: string
    :param fragalysis_dir: the directory containing fragalysis data; within should be target folder, e.g. 'Mpro'
    :type fragalysis_dir: string

    :return: smiles string
    :rtype: string
    """
    dir = os.path.join(fragalysis_dir, target, "aligned")
    fname_part = f"{target}-{fragment}"  # the first part of the filenames/folder names
    protein_path = os.path.join(dir, fname_part, f"{fname_part}_apo-desolv.pdb")

    if return_mol:
        try:
            mol = rdmolfiles.MolFromPDBFile(protein_path)
            return mol
        except OSError as e:
            print("Mol file cannot be found for that fragment.")
            print("Ensure files are saved correctly in Fragalysis format.")
            print(e)
    else:
        return protein_path


def get_files(target, fragment, fragalysis_dir=config_merge.FRAGALYSIS_DATA_DIR):
    """
    Function to get the relevant files for each fragment for filtering.
    File paths are like this: TARGET/aligned/TARGET-FRAGMENT_CHAIN
    e.g. Mpro/aligned/Mpro-x0464_0A

    :param target: the protein target, e.g. 'Mpro'
    :type target: string
    :param fragment: the fragment i.d., e.g. 'x0107'
    :type fragment: string
    :param fragalysis_dir: the directory containing fragalysis data; within should be target folder, e.g. 'Mpro'
    :type fragalysis_dir: string

    :return: mol_file, protein_file
    :rtype: filepaths (strings)
    """
    dir = os.path.join(fragalysis_dir, target, "aligned")
    fname_part = f"{target}-{fragment}"  # the first part of the filenames/folder names
    path = os.path.join(dir, fname_part)
    mol_file = os.path.join(path, f"{fname_part}.mol")
    protein_file = os.path.join(path, f"{fname_part}_apo-desolv.pdb")

    if os.path.exists(mol_file) and os.path.exists(protein_file):
        return mol_file, protein_file
    else:
        print("Cannot find files for that fragment.")
        print("Ensure files are saved correctly in Fragalysis format.")


def get_distance(coord1, coord2):
    """
    Function calculates the distance between two atoms in 3D space.
    Relevant for when two fragments are overlapping.

    :param coord1: 3D atom coordinates
    :type coord1: numpy array
    :param coord2: 3D atom coordinates
    :type coord2: numpy array

    :return: distance between the coordinates
    :rtype: float
    """
    sq = (coord1 - coord2) ** 2
    return np.sqrt(np.sum(sq))


def get_merges(merge_dict):
    """
    Get a list of synthons and matching list of smiles from the merge dictionary
    in the json file. These are then passed into the process_one_smi function
    for filtering.

    :param merge_dict: this is the dictionary containing all the merges
    :type merge_dict: dictionary

    :return: lists of synthons and smiles
    :rtype: lists of strings
    """
    synthons = []
    smiles = []
    for synthon in merge_dict:
        for smi in merge_dict[synthon]:
            if smi:
                synthons.append(synthon)
                smiles.append(smi)
    return synthons, smiles
