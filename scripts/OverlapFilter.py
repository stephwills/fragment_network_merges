"""
Used to filter out compounds that have a large overlap with the protein.
"""

from rdkit import Chem
from rdkit.Chem import rdmolfiles, rdShapeHelpers
import numpy as np
import pandas as pd
from tqdm import tqdm

class OverlapFilter():

    def __init__(self, merges, fragAs, fragBs, proteins, fragsAsFiles=False, proteinAsFile=False, singleProtein=True):
        self.merges = merges  # post-embedding
        self.fragAs = fragAs
        self.fragBs = fragBs
        self.protein = protein
        self.fragsAsFiles = fragsAsFiles  # checks if supplied as mol or file
        self.proteinAsFile = proteinAsFile  # checks if supplied as mol or file
        self.singleProtein = singleProtein  # using a single protein or multiple to calculate distances
        self.distances = []
        self.nonoverlap_mols = []
        self.nonoverlap_indices = []
  
    def merge_protein_distance(self):
        """
        Calculate the distance between the merges and the mols. Accounts for whether a single
        or several proteins are used, as mol(s) already or filenames.
        """
        if self.singleProtein == True:
            if self.proteinAsFile == False:
                for merge in self.merges:
                    distance = rdShapeHelpers.ShapeProtrudeDist(merge, self.protein)
                    self.distances.append(distance)
            else:
                for merge in self.merges:
                    distance = rdShapeHelpers.ShapeProtrudeDist(merge, rdmolfiles.MolFromPDBFile(self.protein))
                    self.distances.append(distance)
        else:
            if self.proteinAsFile == False:
                for merge, protein in zip(self.merges, self.protein):
                    distance = rdShapeHelpers.ShapeProtrudeDist(merge, protein)
                    self.distances.append(distance)
            else:
                for merge, protein in zip(self.merges, self.protein):
                    distance = rdShapeHelpers.ShapeProtrudeDist(merge, rdmolfiles.MolFromPDBFile(protein))
                    self.distances.append(distance)

    def fixed_cutoff_filter(self):
        """
        Rules out molecules that have a >10% overlap with the protein.
        """
        self.fixed_cutoff()
        for i, distance in self.distances:
            if distance >= 0.9:
                self.nonoverlap_indices.append(i)
                self.nonoverlap_mols.append(self.merges[i])
        print(f'Number of molecules: {len(self.nonoverlap_mols)}')
        return self.nonoverlap_indices, self.nonoverlap_mols


class OverlapFilterDf():

    def __init__(self, df, mergeCol='Embedded', fragAsCol='Fragment A file', fragBsCol='Fragment B file', proteinACol='Fragment A protein file', proteinBCol='Fragment B protein file', fragsAsFiles=True, proteinAsFile=True):
        self.df = df
        self.mergeCol = mergeCol
        self.fragAsCol = fragAsCol
        self.fragBsCol = fragBsCol
        self.proteinACol = proteinACol
        self.proteinBCol = proteinBCol
        self.merges = [i for i in df[self.mergeCol]]  # post-embedding
        self.fragAs = [i for i in df[self.fragAsCol]]
        self.fragBs = [i for i in df[self.fragBsCol]]
        self.proteinAs = [i for i in df[self.proteinACol]]
        self.proteinBs = [i for i in df[self.proteinBCol]]
        self.fragsAsFiles = fragsAsFiles  # checks if supplied as mol or file
        self.proteinAsFile = proteinAsFile  # checks if supplied as mol or file
        self.distances = []
        self.nonoverlap_mols = []
        self.nonoverlap_indices = []

    def geometric_mean(self, distA, distB):
        return np.sqrt(distA * distB)

    def merge_protein_distance(self):
        """
        Calculate the distance between the merges and the mols.
        """
        if self.proteinAsFile == False:
            for merge, proteinA, proteinB in tqdm(zip(self.merges, self.proteinAs, self.proteinBs)):
                distanceA = rdShapeHelpers.ShapeProtrudeDist(merge, proteinA)
                distanceB = rdShapeHelpers.ShapeProtrudeDist(merge, proteinB)
                mean = self.geometric_mean(distanceA, distanceB)
                self.distances.append(mean)
        else:
            for merge, proteinA, proteinB in tqdm(zip(self.merges, self.proteinAs, self.proteinBs)):
                distanceA = rdShapeHelpers.ShapeProtrudeDist(merge, rdmolfiles.MolFromPDBFile(proteinA))
                distanceB = rdShapeHelpers.ShapeProtrudeDist(merge, rdmolfiles.MolFromPDBFile(proteinB))
                mean = self.geometric_mean(distanceA, distanceB)
                self.distances.append(mean)

    def fixed_cutoff_filter(self):
        """
        Rules out molecules that have a >10% overlap with the protein.
        """
        self.merge_protein_distance()
        for i, distance in enumerate(self.distances):
            if distance >= 0.9:
                self.nonoverlap_indices.append(i)
                self.nonoverlap_mols.append(self.merges[i])
        print(f'Number of molecules: {len(self.nonoverlap_mols)}')
        new_df = self.df[self.df.index.isin(self.nonoverlap_indices)]
        new_df = new_df.reset_index(drop=True)
        return new_df
