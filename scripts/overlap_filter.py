"""
Used to filter out compounds that have a large overlap with the protein.
"""

from rdkit import Chem
from rdkit.Chem import rdShapeHelpers
import numpy as np

class OverlapFilter():
    """Filters molecules that have >5% overlap with the protein"""

    def __init__(self, merge, proteinA, proteinB):
        self.merge = merge  # molecule post-embedding
        self.proteinA = proteinA  # the protein associated with fragment A
        self.proteinB = proteinB  # the protein associated with fragment B
        self.result = None  # to store the result of the filter

    def geometric_mean(self, distA, distB):
        """
        Calculates the geometric mean between the two distances.

        :param distA: distance between the merge and protein A
        :type distA: float
        :param distB: distance between the merge and protein B
        :type distB: float

        :return: mean distance
        :rtype: float
        """
        return np.sqrt(distA * distB)

    def calc_distances(self):
        """
        Calculate the distance between the merge and the proteins.
        The distance represents the proportion of the volume of the smaller molecule
        that protrudes from the larger molecule (i.e. 1 - overlap).

        :return: distance between the merge and protein A and protein B
        :rtype: float
        """
        distanceA = rdShapeHelpers.ShapeProtrudeDist(self.merge, self.proteinA)
        distanceB = rdShapeHelpers.ShapeProtrudeDist(self.merge, self.proteinB)
        return distanceA, distanceB

    def filter(self):
        """
        Rules out molecules that have a >10% overlap with the protein (i.e. >90%
        protrusion).

        :return: returns 'pass' or 'fail'
        :rtype: string
        """
        distanceA, distanceB = self.calc_distances()  # calculate distances
        mean = self.geometric_mean(distanceA, distanceB)  # calculate mean of distances
        if mean >= 0.9:  # if protrusion > 90% (overlap < 10%)
            self.result = 'pass'
        else:
            self.result = 'fail'
        return self.result
