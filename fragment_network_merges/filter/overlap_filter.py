"""
Used to filter out compounds that clash with the protein.
"""

import argparse
import sys
import os
from typing import Tuple

import numpy as np
from fragment_network_merges.filter.config_filter import config_filter
from fragment_network_merges.filter.generic_filter import Filter_generic
from joblib import Parallel, delayed
from rdkit.Chem import Mol, rdShapeHelpers


class OverlapFilter(Filter_generic):
    def __init__(
        self,
        smis=None,
        synthons=None,
        fragmentA=None,
        fragmentB=None,
        proteinA=None,
        proteinB=None,
        merge=None,
        mols=None,
        names=None,
        work_pair_dir=None,
        out_pair_dir=None,
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
            work_pair_dir,
            out_pair_dir,
        )
        self.results = None

    @staticmethod
    def geometric_mean(distA: float, distB: float) -> float:
        """
        Calculates the geometric mean between the two distances.
        """
        return np.sqrt(distA * distB)

    @staticmethod
    def calc_distances(merge: Mol, proteinA: Mol, proteinB: Mol) -> Tuple[float, float]:
        """
        Calculate the distance between the merge and the proteins. The distance represents the proportion of the volume
        of the smaller molecule that protrudes from the larger molecule (i.e. 1 - overlap).

        :param merge: the merge molecule (after embedding)
        :type merge: rdkit.Chem.Mol
        :param proteinA: protein associated with fragment A
        :type proteinA: rdkit.Chem.Mol
        :param proteinB: protein associated with fragment B
        :type proteinB: rdkit.Chem.Mol

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
        :type merge: rdkit.Chem.Mol
        :param proteinA: protein associated with fragment A
        :type proteinA: rdkit.Chem.Mol
        :param proteinB: protein associated with fragment B
        :type proteinB: rdkit.Chem.Mol
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
        self, cpus: int = config_filter.N_CPUS_FILTER_PAIR, *args
    ) -> Tuple[list, list]:
        """
        Runs the overlap filter on all the SMILES in parallel.

        :param cpus: number of CPUs for parallelization
        :type cpus: int

        :return: list of results (True or False); list of mols (RDKit molecules)
        :rtype: tuple
        """
        self.results = Parallel(n_jobs=cpus, backend="multiprocessing")(
            delayed(self.filter_smi)(mol, self._proteinA, self._proteinB, *args)
            for mol in self.mols
        )
        return self.results, self.mols


def parse_args(args):
    """
    Parse command line arguments
    """
    parser = argparse.ArgumentParser(
        epilog="""
        python filter/overlap_filter.py --input_file data/toFilter.sdf --output_file results.sdf
        --proteinA_file nsp13-x0176_0B_apo-desolv.pdb --proteinB_file nsp13-x0034_0B_apo-desolv.pdb
        --clash_threshold 0.15
        """
    )
    # command line args definitions
    parser.add_argument(
        "-i",
        "--input_file",
        required=True,
        help="input sdf file containing molecules to filter",
    )
    parser.add_argument(
        "-o",
        "--output_file",
        required=True,
        help="output sdf file to write filtered SMILES",
    )
    parser.add_argument(
        "-A", "--proteinA_file", required=True, help="fragment A apo pdb file"
    )
    parser.add_argument(
        "-B", "--proteinB_file", required=True, help="fragment B apo pdb file"
    )
    parser.add_argument(
        "-c",
        "--n_cpus",
        required=False,
        help="number of CPUs",
        type=int,
        default=os.cpu_count(),
    )
    parser.add_argument(
        "-t",
        "--clash_threshold",
        type=float,
        default=0.15,
        help="maximum tolerated clash with protein (proportion of volume of ligand overlapping with protein)",
    )

    return parser.parse_args(args)


def main():
    from fragment_network_merges.filter.generic_squonk import Squonk_generic

    args = parse_args(sys.argv[1:])
    job = Squonk_generic(
        "OverlapFilter",
        args,
        args.input_file,
        args.output_file,
        None,
        None,
        args.proteinA_file,
        args.proteinB_file
    )
    job.execute_job(args.n_cpus, args.clash_threshold)


if __name__ == "__main__":
    main()
