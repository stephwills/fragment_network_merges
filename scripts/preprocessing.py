"""Checks the distance between the fragments"""

import os
import numpy as np
from rdkit import Chem
from rdkit.Chem import rdmolfiles

def get_smiles(target, fragment, chain):
    """
    Function to get the relevant files for each fragment.
    File paths are like this: TARGET/aligned/TARGET-FRAGMENT_CHAIN
    e.g. Mpro/aligned/Mpro-x0464_0A

    :param target: the protein target, e.g. 'Mpro'
    :type target: string
    :param fragment: the fragment i.d., e.g. 'x0107'
    :type fragment: string
    :param chain: the protein chain, e.g. '0A'
    :type chain: string

    :return: smiles string
    :rtype: string
    """
    fname_part = f'{target}-{fragment}_{chain}'  # the first part of the filenames/folder names
    path = os.path.join(target, 'aligned', fname_part)
    smiles_path = os.path.join(path, f'{fname_part}_smiles.txt')
    smiles_file = open(smiles_path)  # open the file to get smiles
    smiles = smiles_file.read()
    smiles_file.close()
    # smiles does not necessarily match what is in the network
    # get canonical smiles by converting to mol and converting back to smiles
    # check this
    mol = Chem.MolFromSmiles(smiles)
    smiles = Chem.MolToSmiles(mol)
    return smiles

def get_mol(target, fragment, chain):
    """
    Function to get the relevant files for each fragment.
    File paths are like this: TARGET/aligned/TARGET-FRAGMENT_CHAIN
    e.g. Mpro/aligned/Mpro-x0464_0A

    :param target: the protein target, e.g. 'Mpro'
    :type target: string
    :param fragment: the fragment i.d., e.g. 'x0107'
    :type fragment: string
    :param chain: the protein chain, e.g. '0A'
    :type chain: string

    :return: smiles string
    :rtype: string
    """
    fname_part = f'{target}-{fragment}_{chain}'  # the first part of the filenames/folder names
    path = os.path.join(target, 'aligned', fname_part)
    mol_file = os.path.join(path, f'{fname_part}.mol')
    mol = rdmolfiles.MolFromMolFile(mol_file)
    return mol

def get_files(target, fragment, chain):
    """
    Function to get the relevant files for each fragment.
    File paths are like this: TARGET/aligned/TARGET-FRAGMENT_CHAIN
    e.g. Mpro/aligned/Mpro-x0464_0A

    :param target: the protein target, e.g. 'Mpro'
    :type target: string
    :param fragment: the fragment i.d., e.g. 'x0107'
    :type fragment: string
    :param chain: the protein chain, e.g. '0A'
    :type chain: string

    :return: smiles_file, mol_file, protein_file
    :rtype: filepaths (strings)
    """
    fname_part = f'{target}-{fragment}_{chain}'  # the first part of the filenames/folder names
    path = os.path.join(target, 'aligned', fname_part)
    smiles_file = os.path.join(path, f'{fname_part}_smiles.txt')
    mol_file = os.path.join(path, f'{fname_part}.mol')
    protein_file = os.path.join(path, f'{fname_part}_apo-desolv.pdb')
    return smiles_file, mol_file, protein_file

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
    returns the average. Used to filter pairs that are too far apart.

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
    
    # get the average of the distances to get an approximation for the dist between the fragments
    avg_distance = sum(distances) / len(distances)
    return avg_distance

def check_fragment_pairs(fragment_pairs, name_pairs, target, chain):
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
        fragmentA = get_mol(target, name_pair[0], chain)
        fragmentB = get_mol(target, name_pair[1], chain)

        # if distance between fragments is >10A, remove this pair
        distance = get_distance_between_fragments(fragmentA, fragmentB)

        if distance < 10:
            filtered_fragment_pairs.append(fragment_pair)
            filtered_name_pairs.append(name_pair)

    return filtered_fragment_pairs, filtered_name_pairs