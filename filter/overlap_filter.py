"""
Used to filter out compounds that have a large overlap with the protein.
"""

from typing import Tuple

import numpy as np
from filter.config_filter import config_filter
from filter.generic_filter import Filter_generic
from joblib import Parallel, delayed
from rdkit.Chem import Mol, rdShapeHelpers


class OverlapFilter(Filter_generic):
    def __init__(
        self,
        smis: list,
        synthons: list,
        fragmentA,
        fragmentB,
        proteinA=None,
        proteinB=None,
        merge=None,
        mols=None,
        names=None,
        pair_working_dir=None,
        pair_output_dir=None,
    ):
        super().__init__(
            smis,
            synthons,
            fragmentA,
            fragmentB,
            proteinA,
            proteinB,
            merge,
            mols,
            names,
            pair_working_dir,
            pair_output_dir,
        )
        self.results = None

    @staticmethod
    def geometric_mean(distA: float, distB: float) -> float:
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

    @staticmethod
    def calc_distances(merge: Mol, proteinA: Mol, proteinB: Mol) -> Tuple[float, float]:
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

    def filter_smi(
        self,
        merge: Mol,
        proteinA: Mol,
        proteinB: Mol,
        clash_dist: float = config_filter.CLASH_DIST,
    ) -> bool:
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
        distanceA, distanceB = self.calc_distances(
            merge, proteinA, proteinB
        )  # calculate distances
        mean = self.geometric_mean(distanceA, distanceB)  # calculate mean of distances
        protrude_dist = 1 - clash_dist
        if mean <= protrude_dist:  # if overlap > clash_dist
            result = False
        else:
            result = True
        return result

    def filter_all(
        self, cpus: int = config_filter.N_CPUS_FILTER_PAIR
    ) -> Tuple[list, list]:
        """
        Runs the descriptor filter on all the SMILES in parallel.

        :param cpus: number of CPUs for parallelization
        :type cpus: int

        :return: list of results (True or False); list of mols (RDKit molecules)
        :rtype: tuple
        """
        self.results = Parallel(n_jobs=cpus, backend="multiprocessing")(
            delayed(self.filter_smi)(mol, self._proteinA, self._proteinB)
            for mol in self.mols
        )
        return self.results, self.mols
