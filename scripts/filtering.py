"""Runs the filtering of the merges for one pair"""

import argparse
import json
import os

import tqdm
from joblib import Parallel, delayed
from rdkit import Chem
from rdkit.Chem import rdmolfiles

from scripts.preprocessing import get_merges, open_json
from scripts.descriptor_filter import descriptor_filter
from scripts.expansion_filter import expansion_filter
from scripts.embedding_filter import embedding_filter
from scripts.overlap_filter import overlap_filter
from scripts.fragmenstein_filter import create_directories, fragmenstein_filter, place_smiles
from scripts.interaction_fp_filter import remove_ligand, ifp_score
from scripts.config import config

def process_one_smi(merge_name, fragmentA, fragmentB, proteinA, proteinB, output_directory, num, smiles, synthon):
    """
    Run all the filters on one smiles. If the smiles fails a filter, the function returns None.
    If it passes all the filters, it returns the smiles.

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

    failed_smiles = None  # stores the smiles if molecule fails (otherwise is None)
    failed_name = None  # stores the name if molecule fails (otherwise is None)
    failed_filter = None  # stores which filter the smiles fails at (if it fails)
    passing_smiles = None  # stores the smiles if molecule is successful (otherwise is None)
    passing_name = None  # stores the name if molecule is successful (otherwise is None)
    if_score = None  # stores ifp score of smiles if successful (otherwise is None)
    result = descriptor_filter(smiles)
    if result == 'fail':
        failed_filter = 'descriptor_filter'
        failed_smiles = smiles
        failed_name = name
        return failed_name, failed_smiles, failed_filter, passing_name, passing_smiles, if_score

    result = expansion_filter(merge_mol, fragmentA_mol, fragmentB_mol, synthon_mol)
    if result == 'fail':
        failed_filter = 'expansion_filter'
        failed_smiles = smiles
        failed_name = name
        return failed_name, failed_smiles, failed_filter, passing_name, passing_smiles, if_score

    embedded, result = embedding_filter(merge_mol, fragmentA_mol, fragmentB_mol, synthon_mol)
    if result == 'fail':
        failed_filter = 'embedding_filter'
        failed_smiles = smiles
        failed_name = name
        return failed_name, failed_smiles, failed_filter, passing_name, passing_smiles, if_score

    result = overlap_filter(embedded, proteinA_mol, proteinB_mol)
    if result == 'fail':
        failed_filter = 'overlap_filter'
        failed_smiles = smiles
        failed_name = name
        return failed_name, failed_smiles, failed_filter, passing_name, passing_smiles, if_score

    try:
        print("launching fragmenstein in %s"%output_directory)
        json_fpath, mol_fpath, pdb_fpath = place_smiles(name, smiles, fragmentA, fragmentB, proteinA, output_directory)
        result = fragmenstein_filter(json_fpath)
        print(result)
        if result == 'fail':
            failed_filter = 'fragmenstein_filter'
            failed_smiles = smiles
            failed_name = name
            return failed_name, failed_smiles, failed_filter, passing_name, passing_smiles, if_score
        else:
            failed_filter = 'no_fail'
            passing_smiles = smiles
            passing_name = name
            pdb_nolig = remove_ligand(pdb_fpath)
            if_score = ifp_score(mol_fpath, fragmentA, fragmentB, pdb_nolig)
            return failed_name, failed_smiles, failed_filter, passing_name, passing_smiles, if_score
    except Exception as e:
            print(e)
            import traceback
            traceback.print_stack()
            failed_filter = 'fragmenstein_filter'
            failed_smiles = smiles
            failed_name = name
            return failed_name, failed_smiles, failed_filter, passing_name, passing_smiles, if_score


def _resolveDataName(fname):
    if os.path.isfile(fname):
        return fname
    else:
        return os.path.join(config.FRAGALYSIS_DATA_DIR, fname)

def main():
    # get file containing merges and all fragment/protein files
    parser = argparse.ArgumentParser()
    parser.add_argument('-f', '--merge_file', help='the json file containing the merges')
    parser.add_argument('-m', '--merge', help='the name of the merge, e.g. x0107_0A_x0434_0A')
    parser.add_argument('-a', '--fragment_A', help='fragment A mol file')
    parser.add_argument('-b', '--fragment_B', help='fragment B mol file')
    parser.add_argument('-p', '--protein_A', help='protein pdb file associated with fragment A')
    parser.add_argument('-q', '--protein_B', help='protein pdb file associated with fragment B')
    parser.add_argument('-o', '--output_directory', help='the directory to write the filtered files to')
    args = parser.parse_args()

    # open json file containing merges
    merges_dict = open_json(args.merge_file)
    synthons, smiles = get_merges(merges_dict)
    print("Number of smiles: %d"%len(smiles))

    # all smiles for each pair are numbered to give unique identifier
    num = range(len(smiles))

    # get args needed for process_one_smi
    fragmentA = _resolveDataName(args.fragment_A)
    fragmentB = _resolveDataName(args.fragment_B)
    proteinA = _resolveDataName(args.protein_A)
    proteinB = _resolveDataName(args.protein_B)
    output_directory = args.output_directory
    merge_name = args.merge

    # create folder to save Fragmenstein files
    create_directories(output_directory)
    # filter all smiles for the pair
    results = Parallel(n_jobs = config.N_CPUS_FILTER_PAIR , backend="multiprocessing")(delayed(process_one_smi)
              (merge_name, fragmentA, fragmentB, proteinA, proteinB, output_directory, n, smi, syn)
              for n, smi, syn in zip(tqdm.tqdm(num), smiles, synthons))

    # if the molecule fails
    failed_name = [res[0] for res in results]
    failed_name = [name for name in failed_name if name]
    failed_smiles = [res[1] for res in results]
    failed_smiles = [smi for smi in failed_smiles if smi]
    failed_filter = [res[2] for res in results]
    failed_filter = [fil for fil in failed_filter if fil != 'no_fail']

    # save info as nested dict
    failed_dict = {}
    for name, smi, fil in zip(failed_name, failed_smiles, failed_filter):
        failed_dict[name] = {'smiles': smi, 'reason': fil}
    
    # save dict of failures to json file
    failures_fname = merge_name + '_failures.json'
    failures_fpath = os.path.join(output_directory, failures_fname)
    with open(failures_fpath, 'w') as f: 
        json.dump(failed_dict, f)

    # if the molecule passes
    passing_name = [res[3] for res in results]
    passing_name = [name for name in passing_name if name]
    passing_smiles = [res[4] for res in results]
    passing_smiles = [smi for smi in passing_smiles if smi]
    if_scores = [res[5] for res in results]
    if_scores = [score for score in if_scores if score]

    # save info as nested dict
    passing_dict = {}
    for name, smi, score in zip(passing_name, passing_smiles, if_scores):
        passing_dict[name] = {'smiles': smi, 'if_score': score}

    # save filtered smiles to json file
    results_fname = merge_name + '_filtered.json'
    results_fpath = os.path.join(output_directory, results_fname)
    with open(results_fpath, 'w') as f:
        json.dump(passing_dict, f)


if __name__ == "__main__":
    main()
