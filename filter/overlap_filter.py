"""
Used to filter out compounds that clash with the protein.
"""

import argparse
import time
from typing import Tuple

import numpy as np
from dm_job_utilities.dm_log import DmLog
from filter.config_filter import config_filter
from filter.generic_filter import Filter_generic
from joblib import Parallel, delayed
from rdkit import Chem
from rdkit.Chem import Mol, rdmolfiles, rdShapeHelpers


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
        self, cpus: int = config_filter.N_CPUS_FILTER_PAIR, **kwargs
    ) -> Tuple[list, list]:
        """
        Runs the descriptor filter on all the SMILES in parallel.

        :param cpus: number of CPUs for parallelization
        :type cpus: int

        :return: list of results (True or False); list of mols (RDKit molecules)
        :rtype: tuple
        """
        self.results = Parallel(n_jobs=cpus, backend="multiprocessing")(
            delayed(self.filter_smi)(mol, self._proteinA, self._proteinB, **kwargs)
            for mol in self.mols
        )
        return self.results, self.mols


def main():
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
        "--clash_threshold",
        type=float,
        default=0.15,
        help="maximum tolerated clash with protein (proportion of volume of ligand overlapping with protein)",
    )

    args = parser.parse_args()

    filter = OverlapFilter()
    DmLog.emit_event("overlap_filter: ", args)

    start = time.time()
    count = 0
    hits = 0
    errors = 0

    proteinA = rdmolfiles.MolFromPDBFile(args.proteinA_file)
    proteinB = rdmolfiles.MolFromPDBFile(args.proteinB_file)

    with Chem.SDWriter(args.output_file) as w:
        with Chem.SDMolSupplier(args.input_file) as suppl:
            for mol in suppl:
                if mol is None:
                    continue
                else:
                    count += 1
                    try:
                        # smi = Chem.MolToSMiles(mol)
                        res = filter.filter_smi(
                            mol, proteinA, proteinB, args.clash_threshold
                        )
                        if res:
                            hits += 1
                            w.write(mol)
                    except Exception as e:
                        DmLog.emit_event(
                            "Failed to process molecule", count, Chem.MolToSmiles(mol)
                        )
                        errors += 1

    end = time.time()
    duration_s = int(end - start)
    if duration_s < 1:
        duration_s = 1

    DmLog.emit_event(
        count, "inputs,", hits, "hits,", errors, "errors.", "Time (s):", duration_s
    )
    DmLog.emit_cost(count)


if __name__ == "__main__":
    main()
