"""
Used for filtering merges according to calculated descriptors.
"""

from typing import Tuple

from filter.config_filter import config_filter
from filter.generic_filter import Filter_generic
from joblib import Parallel, delayed
from rdkit import Chem
from rdkit.Chem import Crippen, rdMolDescriptors


class DescriptorFilter(Filter_generic):
    def __init__(
        self,
        smis: list,
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

    @staticmethod
    def calculate_properties(smi: str) -> int:
        """
        Function to calculate the Lipinski descriptors of the molecule.
        Counts the number of violations of Lipinski's rules.

        :param smi: smiles of the merge
        :type smi: str

        :return: the number of Lipinski rules violated
        :rtype: int
        """
        mol = Chem.MolFromSmiles(smi)

        # calculate the properties
        rotat = rdMolDescriptors.CalcNumRotatableBonds(mol)
        if rotat > 10:  # rule out mols with >10 rb
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
            # print(violations)
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
