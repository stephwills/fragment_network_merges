"""
Used for filtering merges according to calculated descriptors.
"""

import argparse
import time
import sys
from typing import Tuple

from dm_job_utilities.dm_log import DmLog
from filter.config_filter import config_filter
from filter.generic_filter import Filter_generic
from joblib import Parallel, delayed
from rdkit import Chem
from rdkit.Chem import Crippen, rdMolDescriptors


class DescriptorFilter(Filter_generic):
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
    def calculate_properties(
        smi: str,
        rotat_threshold: int = config_filter.ROTAT_THRESHOLD,
        hac_threshold: int = config_filter.HAC_THRESHOLD,
    ) -> int:
        """
        Function to calculate the Lipinski descriptors of the molecule and the number of rotatable bonds.
        If >10 rotatable bonds, will fail filter (violations >1).

        :param smi: smiles of the merge
        :type smi: str
        :param hac_filter: apply hac filter of 15 (for similarity search merges)
        :type hac_filter: bool

        :return: the number of Lipinski rules violated
        :rtype: int
        """
        mol = Chem.MolFromSmiles(smi)

        # calculate the properties
        rotat = rdMolDescriptors.CalcNumRotatableBonds(mol)
        if rotat > rotat_threshold:  # rule out mols with >10 rb
            return 2

        hac = mol.GetNumHeavyAtoms()
        if (
            hac < hac_threshold
        ):  # rule out mols with <15 HAs (for similarity search compounds)
            return 2

        else:
            molw = rdMolDescriptors.CalcExactMolWt(mol)
            alogp = Crippen.MolLogP(mol)
            hba = rdMolDescriptors.CalcNumHBA(mol)
            hbd = rdMolDescriptors.CalcNumHBD(mol)

            # count the number of violations of Lipinski's rules
            violations = 0
            if molw > 500:
                violations += 1
            if alogp > 5:
                violations += 1
            if hba > 10:
                violations += 1
            if hbd > 5:
                violations += 1
            return violations

    def filter_smi(self, smi: str) -> bool:
        """
        Checks whether a SMILES passes or fails Lipinski rules.

        :param smi: smiles of the merge
        :type smi: str

        :return: whether molecule passes (True) or fails (False) filter
        :rtype: bool
        """
        violations = self.calculate_properties(smi)
        if violations <= 1:
            result = True
        else:
            result = False
        return result

    def filter_all(
        self, cpus: int = config_filter.N_CPUS_FILTER_PAIR
    ) -> Tuple[list, None]:
        """
        Runs the descriptor filter on all the SMILES in parallel.

        :param cpus: number of CPUs for parallelization
        :type cpus: int

        :return: list of results (True or False); list of mols (None)
        :rtype: tuple
        """
        self.results = Parallel(n_jobs=cpus, backend="multiprocessing")(
            delayed(self.filter_smi)(smi) for smi in self.smis
        )
        return self.results, self.mols


def parse_args(args):
    """
    Parse command line arguments
    """
    parser = argparse.ArgumentParser(
        epilog="""
    python filter/descriptor_filter.py --input_file data/toFilter.sdf --output_file results.sdf
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
    return parser.parse_args(args)


def main():
    args = parse_args(sys.argv[1:])
    filter = DescriptorFilter()
    DmLog.emit_event("descriptor_filter: ", args)

    start = time.time()
    count = 0
    hits = 0
    errors = 0

    with Chem.SDWriter(args.output_file) as w:
        with Chem.SDMolSupplier(args.input_file) as suppl:
            for mol in suppl:
                if mol is None:
                    continue
                else:
                    count += 1
                    smi = Chem.MolToSmiles(mol)
                    try:
                        res = filter.filter_smi(smi)
                        if res:
                            hits += 1
                            w.write(mol)
                    except Exception as e:
                        DmLog.emit_event("Failed to process molecule", count, smi)
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
