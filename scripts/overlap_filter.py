"""
Used to filter out compounds that have a large overlap with the protein.
"""

from rdkit import Chem
from rdkit.Chem import rdShapeHelpers
import numpy as np


class OverlapFilter():

    def __init__(self, merge, proteinA, proteinB):
        self.merge = merge  # post-embedding
        self.proteinA = proteinA  # the protein associated with fragment A
        self.proteinB = proteinB  # the protein associated with fragment B
        self.result = None  # to store the result of the filter

    def geometric_mean(self, distA, distB):
        """
        Calculates the geometric mean between 
        """
        return np.sqrt(distA * distB)

    def calc_distances(self):
        """
        Calculate the distance between the merge and the proteins.
        """
        distanceA = rdShapeHelpers.ShapeProtrudeDist(self.merge, self.proteinA)
        distanceB = rdShapeHelpers.ShapeProtrudeDist(self.merge, self.proteinB)
        return distanceA, distanceB

    def filter(self):
        """
        Rules out molecules that have a >10% overlap with the protein.
        """
        distanceA, distanceB = self.calc_distances()
        mean = self.geometric_mean(distanceA, distanceB)
        if mean >= 0.9:
            self.result = self.merge
        
        return self.result
