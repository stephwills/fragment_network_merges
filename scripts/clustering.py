"""Cluster fingerprints and create dendrogram"""

import argparse
import os
import pickle
import numpy as np
import matplotlib.pyplot as plt
from rdkit import Chem
from rdkit.Chem import AllChem, DataStructs, rdmolfiles
from scripts.preprocessing import *
from scipy.cluster.hierarchy import dendrogram
from sklearn.cluster import AgglomerativeClustering

def get_smiles_json(json_dir):
    """
    Function to get all the smiles from the directory of jsons containing all the merges

    :param json_dir: directory containing json files
    :type json_dir: directory path
    :return: list of all smiles
    :rtype: list of smiles strings
    """
    fnames = [f for f in os.listdir(json_dir)]  # get all files in directory
    fpaths = [os.path.join(json_dir, fname) for fname in fnames]  # get all file paths
    # get all smiles from the json files
    all_smiles = []
    for fpath in fpaths:
        dict = open_json(fpath)
        _, smiles = get_merges(dict)
        for smi in smiles:
            all_smiles.append(smi)
    return all_smiles

def get_smiles_smi(smi_dir):
    """
    Function to get all the smiles from the directory of jsons containing all the filtered merges

    :param json_dir: directory containing json files
    :type json_dir: directory path
    :return: list of all smiles
    :rtype: list of smiles strings
    """
    fnames = [f for f in os.listdir(smi_dir)]  # get all files in directory
    fpaths = [os.path.join(smi_dir, fname) for fname in fnames]  # get all file paths
    # get all smiles from the json files
    all_smiles = []
    for fpath in fpaths:
        filtered = open_json(fpath)
        if len(filtered) != 0:
            for smi in filtered:
                all_smiles.append(smi)
    return all_smiles

def get_fps(smiles):
    """
    Function to convert all smiles to Morgan fingerprints (radius; 1024 bits)

    :param smiles: all smiles to convert
    :type smiles: list of smiles strings
    :return: list of fingerprints
    :rtype: list of bit vectors
    """
    # convert smiles to mol and convert to fp
    fps = [AllChem.GetMorganFingerprintAsBitVect(Chem.MolFromSmiles(smi), 2, nBits=1024) for smi in smiles]
    return fps

def tanimoto(x, y):
    """
    Calculate the Tanimoto similarity between two fingerprints (for distance matrix)

    :param x: Morgan fingerprint
    :type x: bit vector
    :param y: Morgan fingerprint
    :type y: bit vector
    :return: Tanimoto metric
    :rtype: float
    """
    return DataStructs.cDataStructs.TanimotoSimilarity(x, y)

def get_dist_matrix(fps):
    """
    Function calculates the distance matrix for all fingerprints in the compound set for clustering

    :param fps: list of fingerprints
    :type fps: list of bit vectors
    :return: distance matrix
    :rtype: 2D numpy array
    """
    nfps = len(fps)
    dists = np.empty([nfps, nfps])  # create empty 2D numpy array to store distances
    for i, fp1 in enumerate(fps):
        for j, fp2 in enumerate(fps):
            sim = tanimoto(fp1, fp2)  # calculate similarity (Tanimoto)
            dists[i, j] = 1 - sim  # distance is 1 - similarity
    return dists

def cluster_fps(dists):
    """
    Function clusters the distance matrix and returns the model

    :param dists: distance matrix
    :type dists: 2D numpy array
    :return: fitted model
    :rtype: sklearn model
    """
    # complete linkage uses the maximum distances between all observations of the two sets
    # affinity set as 'precomputed' to indicate distance matrix used
    clustering = AgglomerativeClustering(affinity='precomputed', distance_threshold=0, n_clusters=None, linkage='complete')
    clustering = clustering.fit(dists)
    return clustering

def plot_dendrogram(model, **kwargs):
    """
    Creates linkage matrix and plots dendrogram

    :param model: clustering model
    :type model: sklearn fitted model
    """
    # create the counts of samples under each node
    counts = np.zeros(model.children_.shape[0])
    n_samples = len(model.labels_)
    for i, merge in enumerate(model.children_):
        current_count = 0
        for child_idx in merge:
            if child_idx < n_samples:
                current_count += 1  # leaf node
            else:
                current_count += counts[child_idx - n_samples]
        counts[i] = current_count

    linkage_matrix = np.column_stack([model.children_, model.distances_,
                                      counts]).astype(float)

    # plot the corresponding dendrogram
    # truncate_mode: used to condense the dendrogram
    # level means no more than p levels are used
    dendrogram(linkage_matrix, **kwargs)

def run_clustering():
    # get arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('-d', '--json_directory', help='directory containing the jsons with merges')
    parser.add_argument('-s', '--smiles_directory', help='directory containing the jsons with the filtered compounds')
    parser.add_argument('-o', '--output_directory', help='output directory to save files to')
    parser.add_argument('-f', '--output_file', help='what to name the files (without file extension)')
    args = parser.parse_args()

    # execute functions
    if args.json_directory:
        json_dir = args.json_directory
        smiles = get_smiles_json(json_dir)
    if args.smiles_directory:
        smi_dir = args.smiles_directory
        smiles = get_smiles_smi(smi_dir)
    output_dir = args.output_directory
    output_file = args.output_file

    fps = get_fps(smiles)
    dists = get_dist_matrix(fps)
    model = cluster_fps(dists)

    # plot model
    plt.figure(figsize=(14,10))
    plt.title('Hierarchical clustering dendrogram', fontsize=25)
    plot_dendrogram(model, truncate_mode='level', p=3)  # no more than 5 levels of tree displayed
    plt.xlabel('Number of molecules in node', fontsize=20)
    plt.xticks(fontsize=15)
    plt.yticks(fontsize=15)
    fname_save = output_file + '.png'
    file_to_save = os.path.join(output_dir, fname_save)
    plt.savefig(file_to_save)

    # save model
    mname_save = output_file + '.pkl'
    model_to_save = os.path.join(output_dir, mname_save)
    pickle.dump(model, open(model_to_save, "wb"))

if __name__ == "__main__":
    run_clustering()
