"""
Used to filter out compounds that look like elaborations rather than merges.
"""

from typing import Tuple

from filter.config_filter import config_filter
from filter.generic_filter import Filter_generic
from joblib import Parallel, delayed
from rdkit import Chem
from rdkit.Chem import AllChem, Mol, rdFMCS

# Function has been massively simplified; may need more testing


class ExpansionFilter(Filter_generic):
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
    def _calc_mcs(fragment: Mol, merge: Mol) -> Mol:
        """Function checks if there is a MCS between fragment A, B and the merge.
        Returns the molecule if there is and returns None if not.

        :param fragment:
        :param merge:

        :return: the MCS molecule or None
        :rtype: RDKit molecule or None
        """
        mcs = rdFMCS.FindMCS([fragment, merge], completeRingsOnly=True)
        print('MCS', mcs)
        mcs_mol = Chem.MolFromSmarts(mcs.smartsString)
        Chem.SanitizeMol(mcs_mol)
        return mcs_mol

    def filter_smi(
        self,
        smiles: str,
        fragmentA: Mol,
        fragmentB: Mol,
        volume_threshold: float = config_filter.FRAGA_VOLUME,
    ) -> bool:
        """
        Checks if molecule resemble expansions of fragment A (where fragment B is not contributing
        anything unique to the merge). Calculates the MCS between fragment A and the merge. If the MCS
        accounts for >90% the volume of the merge, the function checks whether the remaining part
        is from fragment B. If not, the molecule is rejected. If it is from fragment B, the molecule
        passes.

        :param smiles: smiles string of the merge
        :type smiles: str
        :param fragmentA: fragment A molecule
        :type fragmentA: RDKit molecule
        :param fragmentB: fragment B molecule
        :type fragmentB: RDKit molecule
        :volume_threshold: the max proportion of volume accounted for by fragment A
        :type volume_threshold: float

        :return: whether molecule passes (True) or fails (False) filter
        :rtype: bool
        """
        merge = Chem.MolFromSmiles(smiles)
        mcs = self._calc_mcs(fragmentA, merge)

        AllChem.EmbedMolecule(merge)
        AllChem.EmbedMolecule(mcs)
        merge_vol = AllChem.ComputeMolVolume(merge)
        mcs_vol = AllChem.ComputeMolVolume(mcs)
        print((mcs_vol / merge_vol))
        if (mcs_vol / merge_vol) > volume_threshold:  # if mcs account for > 90%
            print('Volume > 90%')
            non_mcs = AllChem.DeleteSubstructs(merge, mcs)
            match = fragmentB.HasSubstructMatch(non_mcs)
            print(match)
            if match:
                result = True
            else:
                result = False

        else:
            print('Volume < 90%')
            result = True

        return result

    def filter_all(
        self, cpus: int = config_filter.N_CPUS_FILTER_PAIR
    ) -> Tuple[list, None]:
        """
        Runs the expansion filter on all the SMILES in parallel.

        :param cpus: number of CPUs for parallelization
        :type cpus: int

        :return: list of results (True or False); list of mols (None)
        :rtype: tuple
        """
        self.results = Parallel(n_jobs=cpus, backend="multiprocessing")(
            delayed(self.filter_smi)(smi, self._fragmentA, self._fragmentB)
            for smi in self.smis
        )
        return self.results, self.mols
