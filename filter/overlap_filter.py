"""
Used to filter out compounds that have a large overlap with the protein.
"""

import numpy as np

from joblib import Parallel, delayed
from rdkit.Chem import rdShapeHelpers
from typing import Tuple

from filter.config_filter import config_filter
from filter.generic_filter import Filter_generic


class OverlapFilter(Filter_generic):

    def __init__(self, smis: list, synthons: list, fragmentA, fragmentB, proteinA, proteinB, merge, mols, names):
        super().__init__(smis, synthons, fragmentA, fragmentB, proteinA, proteinB, merge, mols, names)
        self.results = None

    def geometric_mean(self, distA: float, distB: float) -> float:
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

    def calc_distances(self, merge, proteinA, proteinB) -> Tuple[float, float]:
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

    def filter_smi(self, merge, proteinA, proteinB, clash_dist: float = config_filter.CLASH_DIST) -> bool:
        """
        Rules out molecules that have a >10% overlap with the protein (i.e. >90%
        protrusion).

        :param merge: the merge molecule (after embedding)
        :type merge: RDKit molecule
        :param proteinA: protein associated with fragment A
        :type proteinA: RDKit molecule
        :param proteinB: protein associated with fragment B
        :type proteinB: RDKit molecule
        :param clash_dist: the threshold for the amount of overlap allowed with the protein
        :type clash_dist: float

        :return: whether molecule passes (True) or fails (False) filter
        :rtype: bool
        """
        distanceA, distanceB = self.calc_distances(merge, proteinA, proteinB)  # calculate distances
        mean = self.geometric_mean(distanceA, distanceB)  # calculate mean of distances
        protrude_dist = 1 - clash_dist
        if mean <= protrude_dist:  # if overlap > clash_dist
            result = False
        else:
            result = True
        return result

    def filter_all(self, cpus: int = config_filter.N_CPUS_FILTER_PAIR) -> Tuple[list, list]:
        """
        Runs the descriptor filter on all the SMILES in parallel.

        :param cpus: number of CPUs for parallelization
        :type cpus: int

        :return: list of results (True or False); list of mols (RDKit molecules)
        :rtype: tuple
        """
        self.results = Parallel(n_jobs=cpus, backend='multiprocessing') \
            (delayed(self.filter_smi)(mol, self._proteinA, self._proteinB) for mol in self.mols)
        return self.results, self.mols
