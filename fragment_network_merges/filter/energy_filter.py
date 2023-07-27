"""
Used for filtering out poses with bad energies.
"""

import argparse
import os
import sys
from typing import Tuple

from fragment_network_merges.filter.config_filter import config_filter
from fragment_network_merges.filter.generic_filter import Filter_generic
from joblib import Parallel, delayed
from rdkit import RDLogger
from rdkit.Chem import Mol
from fragment_network_merges.utils.filter_utils import calc_energy, calc_unconstrained_energy

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
        try:
            const_energy = calc_energy(merge)
            print("constrained", const_energy)
            # energy of avg unconstrained conformation
            unconst_energy = calc_unconstrained_energy(merge, n_conf)
            print("unconstrained", unconst_energy)
            # if the energy of the constrained conformation is less, then pass filter
            if (const_energy / unconst_energy) >= energy_threshold:
                result = False
            else:
                result = True
            return result
        except Exception as e:
            print(e)
            return False

    def filter_all(
        self, cpus: int = config_filter.N_CPUS_FILTER_PAIR, *args
    ) -> Tuple[list, list]:
        """
        Runs the energy filter on all the SMILES in parallel.

        :param cpus: number of CPUs for parallelization
        :type cpus: int

        :return: list of results (True or False); list of mols (RDKit molecules)
        :rtype: tuple
        """
        self.results = Parallel(n_jobs=cpus, backend="multiprocessing")(
            delayed(self.filter_smi)(mol, *args) for mol in self.mols
        )
        return self.results, self.mols


def parse_args(args):
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
        "-c",
        "--n_cpus",
        required=False,
        help="number of CPUs",
        type=int,
        default=os.cpu_count(),
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

    return parser.parse_args(args)


def main():
    from fragment_network_merges.filter.generic_squonk import Squonk_generic

    args = parse_args(sys.argv[1:])
    job = Squonk_generic("EnergyFilter", args, args.input_file, args.output_file)
    job.execute_job(args.n_cpus, args.energy_threshold, args.n_conformations)


if __name__ == "__main__":
    main()
