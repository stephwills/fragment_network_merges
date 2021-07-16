"""
Used to filter out compounds that have a large overlap with the protein.
"""

import numpy as np

from rdkit import Chem
from rdkit.Chem import rdShapeHelpers

def geometric_mean(distA, distB):
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

def calc_distances(merge, proteinA, proteinB):
    """
    Calculate the distance between the merge and the proteins.
    The distance represents the proportion of the volume of the smaller molecule
    that protrudes from the larger molecule (i.e. 1 - overlap).

    :param merge: the merge molecule (after embedding)
    :type merge: RDKit molecule
    :param proteinA: protein associated with fragment A
    :type proteinA: RDKit mol from pdb file
    :param proteinB: protein associated with fragment B
    :type proteinB: RDKit mol from pdb file

    :return: distance between the merge and protein A and protein B
    :rtype: float
    """
    distanceA = rdShapeHelpers.ShapeProtrudeDist(merge, proteinA)
    distanceB = rdShapeHelpers.ShapeProtrudeDist(merge, proteinB)
    return distanceA, distanceB

def overlap_filter(merge, proteinA, proteinB):
    """
    Rules out molecules that have a >10% overlap with the protein (i.e. >90%
    protrusion).

    :param merge: the merge molecule (after embedding)
    :type merge: RDKit molecule
    :param proteinA: protein associated with fragment A
    :type proteinA: RDKit mol from pdb file
    :param proteinB: protein associated with fragment B
    :type proteinB: RDKit mol from pdb file

    :return: returns 'pass' or 'fail'
    :rtype: string
    """
    distanceA, distanceB = calc_distances(merge, proteinA, proteinB)  # calculate distances
    mean = geometric_mean(distanceA, distanceB)  # calculate mean of distances
    if mean >= 0.85:  # if protrusion > 85% (overlap < 15%)
        result = 'pass'
    else:
        result = 'fail'
    return result
