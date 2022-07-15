"""
Used for filtering out poses with bad energies.
"""

import argparse
import time
from typing import Tuple

from dm_job_utilities.dm_log import DmLog
from filter.config_filter import config_filter
from filter.generic_filter import Filter_generic
from joblib import Parallel, delayed
from rdkit import Chem, RDLogger
from rdkit.Chem import Mol
from utils.filter_utils import calc_energy, calc_unconstrained_energy


RDLogger.DisableLog("rdApp.*")

class EnergyFilter(Filter_generic):
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

    def filter_smi(
        self,
        merge: Mol,
        energy_threshold: float = config_filter.ENERGY_THRESHOLD,
        n_conf: int = config_filter.N_CONFORMATIONS,
    ) -> Tuple[bool, None]:
        """
        Runs the energy filter. Calculates the ratio between the energy of the constrained conformation and the
        average for multiple unconstrained conformations, and rules out those above a certain threshold.

        :param merge: merge molecule
        :type merge: rdkit.Chem.Mol
        :param energy_threshold: threshold of constrained energy versus unconstrained energy ratio
        :type energy_threshold: float
        :param n_conf: number of unconstrained conformations to generate
        :type n_conf: int

        :return: whether molecule passes (True) or fails (False) filter
        :rtype: bool
        """
        const_energy = calc_energy(merge)
        # energy of avg unconstrained conformation
        unconst_energy = calc_unconstrained_energy(merge, n_conf)
        # if the energy of the constrained conformation is less, then pass filter
        if (const_energy / unconst_energy) >= energy_threshold:
            result = False
        else:
            result = True
        return result

    def filter_all(
        self, cpus: int = config_filter.N_CPUS_FILTER_PAIR, **kwargs
    ) -> Tuple[list, list]:
        """
        Runs the energy filter on all the SMILES in parallel.

        :param cpus: number of CPUs for parallelization
        :type cpus: int

        :return: list of results (True or False); list of mols (RDKit molecules)
        :rtype: tuple
        """
        self.results = Parallel(n_jobs=cpus, backend="multiprocessing")(
            delayed(self.filter_smi)(mol, **kwargs) for mol in self.mols
        )
        return self.results, self.mols


def main():
    parser = argparse.ArgumentParser(
        epilog="""
    python filter/energy_filter.py --input_file data/toFilter.sdf --output_file results.sdf
    --energy_threshold 10 --n_conformations 50
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
        "-e",
        "--energy_threshold",
        type=float,
        default=10.0,
        help="threshold for fold difference in energy of unconstrained versus constrained conformations",
    )
    parser.add_argument(
        "-n",
        "--n_conformations",
        type=int,
        default=50,
        help="number of unconstrained conformations to generate",
    )

    args = parser.parse_args()

    filter = EnergyFilter()
    DmLog.emit_event("energy_filter: ", args)

    start = time.time()
    count = 0
    hits = 0
    errors = 0
    print(args.energy_threshold)
    print(type(args.energy_threshold))
    with Chem.SDWriter(args.output_file) as w:
        with Chem.SDMolSupplier(args.input_file) as suppl:
            for mol in suppl:
                if mol is None:
                    continue
                else:
                    print(mol)
                    count += 1
                    try:

                        res = filter.filter_smi(
                            mol,
                            args.energy_threshold,
                            args.n_conformations,
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
