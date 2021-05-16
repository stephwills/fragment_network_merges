"""Functions for preprocessing files/fragments before running merge generation/filtering"""

import os
import json
import numpy as np
from rdkit import Chem
from rdkit.Chem import rdmolfiles

def open_json(file):
    """
    Function to open json file
    """
    f = open(file, 'r')
    data = json.load(f)
    return data

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

def get_smiles(target, fragment):
    """
    Function to get the SMILES for each fragment.
    File paths are like this: TARGET/aligned/TARGET-FRAGMENT_CHAIN
    e.g. Mpro/aligned/Mpro-x0464_0A

    :param target: the protein target, e.g. 'Mpro'
    :type target: string
    :param fragment: the fragment i.d., e.g. 'x0107_0A'
    :type fragment: string

    :return: smiles string
    :rtype: string
    """
    fname_part = f'{target}-{fragment}'  # the first part of the filenames/folder names
    path = os.path.join(target, 'aligned', fname_part)
    smiles_path = os.path.join(path, f'{fname_part}_smiles.txt')
    smiles_file = open(smiles_path)  # open the file to get smiles
    smiles = smiles_file.read()
    smiles_file.close()
    # smiles does not necessarily match what is in the network
    # get canonical smiles by converting to mol and converting back to smiles
    mol = Chem.MolFromSmiles(smiles)
    smiles = Chem.MolToSmiles(mol)
    return smiles

def get_mol(target, fragment):
    """
    Function to get the mol for each fragment.
    File paths are like this: TARGET/aligned/TARGET-FRAGMENT_CHAIN
    e.g. Mpro/aligned/Mpro-x0464_0A

    :param target: the protein target, e.g. 'Mpro'
    :type target: string
    :param fragment: the fragment i.d., e.g. 'x0107'
    :type fragment: string

    :return: smiles string
    :rtype: string
    """
    fname_part = f'{target}-{fragment}'  # the first part of the filenames/folder names
    path = os.path.join(target, 'aligned', fname_part)
    mol_file = os.path.join(path, f'{fname_part}.mol')
    mol = rdmolfiles.MolFromMolFile(mol_file)
    return mol

def get_files(target, fragment):
    """
    Function to get the relevant files for each fragment for filtering.
    File paths are like this: TARGET/aligned/TARGET-FRAGMENT_CHAIN
    e.g. Mpro/aligned/Mpro-x0464_0A

    :param target: the protein target, e.g. 'Mpro'
    :type target: string
    :param fragment: the fragment i.d., e.g. 'x0107'
    :type fragment: string

    :return: mol_file, protein_file
    :rtype: filepaths (strings)
    """
    fname_part = f'{target}-{fragment}'  # the first part of the filenames/folder names
    path = os.path.join(target, 'aligned', fname_part)
    mol_file = os.path.join(path, f'{fname_part}.mol')
    protein_file = os.path.join(path, f'{fname_part}_apo-desolv.pdb')
    return mol_file, protein_file

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

def get_distance_between_fragments(fragmentA, fragmentB):
    """
    Function to calculate the distance between two fragments. Calculates
    the distance between each pair of atoms in the two fragments, and
    returns the shortest distance. Used to filter pairs that are too far apart.

    :param fragmentA: fragment A molecule
    :type fragmentA: RDKit molecule
    :param fragmentB: fragment B molecule
    :type fragmentB: RDKit molecule

    :return: the average distance between the two fragments
    :rtype: float
    """
    # get the conformers from the fragments
    confA = fragmentA.GetConformer()
    confB = fragmentB.GetConformer()

    # store all the distances between each pair of atoms in the two fragments
    distances = []

    # for each atom in A, calculate the distance with each atom in B
    # append to distances
    for i in range(confA.GetNumAtoms()):
        posA = np.array(confA.GetAtomPosition(i))  # get 3D coordinates
        for j in range(confB.GetNumAtoms()):
            posB = np.array(confB.GetAtomPosition(j))  # get 3D coordinates
            distance = get_distance(posA, posB)  # calculate distance
            distances.append(distance)  # append to list

    # get the shortest distance between the fragments
    return min(distances)

def check_fragment_pairs(fragment_pairs, name_pairs, target):
    """
    Function to filter the list of fragment pairs for those that are close
    enough to merge.

    :param fragment_pairs: list of tuples of fragment pairs
    :type fragment_pairs: list of tuples
    :param name_pairs: list of tuples of fragment names
    :type name_pairs: list of tuples

    :return: the two lists filtered
    :rtype: two list of tuples
    """
    filtered_fragment_pairs = []
    filtered_name_pairs = []
    for fragment_pair, name_pair in zip(fragment_pairs, name_pairs):
        fragmentA = get_mol(target, name_pair[0])
        fragmentB = get_mol(target, name_pair[1])

        # if distance between fragments is >5A, remove this pair
        distance = get_distance_between_fragments(fragmentA, fragmentB)

        if distance < 5:
            filtered_fragment_pairs.append(fragment_pair)
            filtered_name_pairs.append(name_pair)

    # write fragment pairs list to json file
    filename = os.path.join('data', 'fragment_pairs.json')
    with open(filename, 'w') as f:
        json.dump(filtered_name_pairs, f)

    return filtered_fragment_pairs, filtered_name_pairs
