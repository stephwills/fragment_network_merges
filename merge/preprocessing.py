"""Functions for preprocessing files/fragments before running merge generation/filtering"""

import os
import json
import numpy as np

from rdkit import Chem
from rdkit.Chem import rdmolfiles
from merge.config_merge import config_merge


def load_json(fname):
    """
    Function to open json file
    """
    with open(fname, 'r') as f:
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
    dir = os.path.join(fragalysis_dir, target, 'aligned')
    fname_part = f'{target}-{fragment}'  # the first part of the filenames/folder names
    path = os.path.join(dir, fname_part)
    smiles_path = os.path.join(path, f'{fname_part}_smiles.txt')

    try:
        with open(smiles_path) as smiles_file:  # open the file to get smiles
            smiles = smiles_file.read()
    except OSError as e:
        print('SMILES file cannot be found for that fragment.')
        print('Ensure files are saved correctly in Fragalysis format.')
        print(e)

    # smiles does not necessarily match what is in the network
    # get canonical smiles by converting to mol and converting back to smiles
    mol = Chem.MolFromSmiles(smiles)
    smiles = Chem.MolToSmiles(mol)
    return smiles


def get_mol(target, fragment, return_mol=False, fragalysis_dir=config_merge.FRAGALYSIS_DATA_DIR):
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
    :rtype: string
    """
    dir = os.path.join(fragalysis_dir, target, 'aligned')
    fname_part = f'{target}-{fragment}'  # the first part of the filenames/folder names
    path = os.path.join(dir, fname_part)
    mol_path = os.path.join(path, f'{fname_part}.mol')

    if return_mol:
        try:
            mol = rdmolfiles.MolFromMolFile(mol_path)
            return mol
        except OSError as e:
            print('Mol file cannot be found for that fragment.')
            print('Ensure files are saved correctly in Fragalysis format.')
            print(e)
    else:
        return mol_path


def get_protein(target, fragment, fragalysis_dir=config_merge.FRAGALYSIS_DATA_DIR, return_mol=False):
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
    dir = os.path.join(fragalysis_dir, target, 'aligned')
    fname_part = f'{target}-{fragment}'  # the first part of the filenames/folder names
    protein_path = os.path.join(dir, fname_part, f'{fname_part}_apo-desolv.pdb')

    if return_mol:
        try:
            mol = rdmolfiles.MolFromPDBFile(protein_path)
            return mol
        except OSError as e:
            print('Mol file cannot be found for that fragment.')
            print('Ensure files are saved correctly in Fragalysis format.')
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
    dir = os.path.join(fragalysis_dir, target, 'aligned')
    fname_part = f'{target}-{fragment}'  # the first part of the filenames/folder names
    path = os.path.join(dir, fname_part)
    mol_file = os.path.join(path, f'{fname_part}.mol')
    protein_file = os.path.join(path, f'{fname_part}_apo-desolv.pdb')

    if os.path.exists(mol_file) and os.path.exists(protein_file):
        return mol_file, protein_file
    else:
        print('Cannot find files for that fragment.')
        print('Ensure files are saved correctly in Fragalysis format.')


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


def check_fragment_pairs(fragment_pairs, name_pairs, target, max_dist=config_merge.MAX_FRAG_DIST,
                         working_dir=config_merge.WORKING_DIR):
    """
    Function to filter the list of fragment pairs for those that are close
    enough to merge.

    :param fragment_pairs: list of tuples of fragment pairs
    :type fragment_pairs: list of tuples
    :param name_pairs: list of tuples of fragment names
    :type name_pairs: list of tuples
    :param target: the protein target, e.g. 'Mpro'
    :type target: string
    :param max_dist: maximum distance between fragments to be considered for merging
    :type max_dist: float
    :param working_dir: the current working directory to save intermediate files
    :type working_dir: string

    :return: the two lists filtered
    :rtype: two list of tuples
    """
    filtered_fragment_pairs = []
    filtered_name_pairs = []

    # filter the fragment pairs by distance between them
    for fragment_pair, name_pair in zip(fragment_pairs, name_pairs):
        fragmentA = get_mol(target, name_pair[0], True)
        fragmentB = get_mol(target, name_pair[1], True)

        # if distance between fragments is >max_dist, remove this pair
        distance = get_distance_between_fragments(fragmentA, fragmentB)

        if distance <= max_dist:
            filtered_fragment_pairs.append(fragment_pair)
            filtered_name_pairs.append(name_pair)

    # write fragment pairs list to json file
    if not os.path.exists(working_dir):
        os.mkdir(working_dir)
    filename = os.path.join(working_dir, f'{target}_pairs.json')

    # check if fragment pairs list exists for this target in the working directory
    if os.path.exists(filename):
        existing_name_pairs = load_json(filename)
        if sorted(existing_name_pairs) == sorted(filtered_name_pairs):  # check if fragment pairs are equivalent
            print(f'Fragment pairs have already been enumerated for this set of fragments against {target}')
            return filtered_fragment_pairs, filtered_name_pairs
        else:
            print(f'A different set of fragment pairs have already been enumerated against {target}')

    with open(filename, 'w') as f:
        json.dump(filtered_name_pairs, f)

    return filtered_fragment_pairs, filtered_name_pairs


def get_pair_dict(fragment_pairs):
    """
    Create dictionary with fragment A as key and all the fragment Bs used for expansion as the vals (in list).
    """
    pair_dict = {}
    for pair in fragment_pairs:
        fA, fB = pair[0], pair[1]
        if fA in pair_dict:
            pair_dict[fA].append(fB)
        else:
            pair_dict[fA] = [fB]
    return pair_dict


def check_merges_run(smiles_pairs, name_pairs, output_dir):
    """
    Checks if a json file already exists in the output directory
    for the file (avoid re-running queries).
    """
    if output_dir:
        if not os.path.exists(output_dir):
            os.mkdir(output_dir)

    already_run = []  # record pairs already run
    for i, pair in enumerate(name_pairs):
        merge = pair[0] + '_' + pair[1]
        if output_dir:
            filename = merge + '.json'
            filepath = os.path.join(output_dir, filename)
            if os.path.isfile(filepath):  # if file exists for merge pair then remove from list
                smiles_pairs.pop(i)
                name_pairs.pop(i)
                already_run.append(merge)
                print('The following merges have already been run', already_run)
                print(f'{len(name_pairs)} merge pairs remaining')

        return smiles_pairs, name_pairs
