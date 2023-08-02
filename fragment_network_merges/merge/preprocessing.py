"""
Functions for preprocessing before running merge generation
"""

import json
import os

import numpy as np
from fragment_network_merges.merge.config_merge import config_merge
from rdkit import Chem
from rdkit.Chem import rdFMCS
from fragment_network_merges.utils.utils import get_distance, get_mol, load_json


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


def check_fragment_pairs(
    fragment_pairs,
    name_pairs,
    target,
    max_dist=config_merge.MAX_FRAG_DIST,
    working_dir=config_merge.WORKING_DIR,
):
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
    filename = os.path.join(working_dir, f"{target}_pairs.json")

    # check if fragment pairs list exists for this target in the working directory
    if os.path.exists(filename):
        existing_name_pairs = load_json(filename)
        if sorted(existing_name_pairs) == sorted(
            filtered_name_pairs
        ):  # check if fragment pairs are equivalent
            print(
                f"Fragment pairs have already been enumerated for this set of fragments against {target}"
            )
            return filtered_fragment_pairs, filtered_name_pairs
        else:
            print(
                f"A different set of fragment pairs have already been enumerated against {target}"
            )

    with open(filename, "w") as f:
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
        merge = pair[0] + "_" + pair[1]
        if output_dir:
            filename = merge + ".json"
            filepath = os.path.join(output_dir, filename)
            if os.path.isfile(
                filepath
            ):  # if file exists for merge pair then remove from list
                smiles_pairs.pop(i)
                name_pairs.pop(i)
                already_run.append(merge)
                print("The following merges have already been run", already_run)
                print(f"{len(name_pairs)} merge pairs remaining")

    return smiles_pairs, name_pairs


def calculate_msd(molA, molB, mapping):
    confA = molA.GetConformer()
    confB = molB.GetConformer()
    return sum([(confA.GetAtomPosition(a).x - confB.GetAtomPosition(b).x) ** 2 +
                (confA.GetAtomPosition(a).y - confB.GetAtomPosition(b).y) ** 2 +
                (confA.GetAtomPosition(a).z - confB.GetAtomPosition(b).z) ** 2 for a, b in mapping])


def calculate_rmsd(molA, molB, mapping) -> float:
    return (calculate_msd(molA, molB, mapping) / len(mapping)) ** 0.5


def check_similarity(mol1, mol2, num_atoms_not_mcs=config_merge.NUM_ATOMS_NOT_MCS, rmsd_threshold=config_merge.MCS_RMSD_THRESHOLD):
    mcs = Chem.MolFromSmarts(rdFMCS.FindMCS([mol1, mol2], completeRingsOnly=True).smartsString)
    mcs_atoms = mcs.GetNumAtoms()
    mol1_atoms = mol1.GetNumAtoms() - mcs_atoms
    mol2_atoms = mol2.GetNumAtoms() - mcs_atoms
    if mol1_atoms > num_atoms_not_mcs or mol2_atoms > num_atoms_not_mcs:
        return True

    mapping = []
    mol1_matches = mol1.GetSubstructMatch(mcs)
    mol2_matches = mol2.GetSubstructMatch(mcs)
    for mol1_match, mol2_match in zip(mol1_matches, mol2_matches):
        mapping.append((mol1_match, mol2_match))

    rmsd = calculate_rmsd(mol1, mol2, mapping)
    if rmsd <= rmsd_threshold:
        return False

    else:
        return True


def check_too_similar(fragment_pairs, name_pairs, target, num_atoms_not_mcs=config_merge.NUM_ATOMS_NOT_MCS,
                      rmsd_threshold=config_merge.MCS_RMSD_THRESHOLD):
    filtered_fragment_pairs = []
    filtered_name_pairs = []

    for fragment_pair, name_pair in zip(fragment_pairs, name_pairs):
        fragmentA = get_mol(target, name_pair[0], True)
        fragmentB = get_mol(target, name_pair[1], True)

        check_similar = check_similarity(fragmentA, fragmentB, num_atoms_not_mcs, rmsd_threshold)
        if check_similar:
            filtered_fragment_pairs.append(fragment_pair)
            filtered_name_pairs.append(name_pair)

    return filtered_fragment_pairs, filtered_name_pairs
