"""
Used to filter out compounds that look like elaborations rather than merges.
"""

from typing import Tuple

from filter.config_filter import config_filter
from filter.embedding_filter import remove_xe
from filter.generic_filter import Filter_generic
from joblib import Parallel, delayed
from rdkit import Chem
from rdkit.Chem import Lipinski, Mol, rdFMCS


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
    def _sanitization(mol: Mol):
        """
        Do a partial sanitization of mol if normal sanitization not possible (e.g. due to valence errors after atom
        removal or for MCS. Code from RDKit cookbook.
        """
        try:
            Chem.SanitizeMol(mol)

        except:
            mol.UpdatePropertyCache(strict=False)
            Chem.SanitizeMol(
                mol,
                Chem.SanitizeFlags.SANITIZE_FINDRADICALS
                | Chem.SanitizeFlags.SANITIZE_KEKULIZE
                | Chem.SanitizeFlags.SANITIZE_SETAROMATICITY
                | Chem.SanitizeFlags.SANITIZE_SETCONJUGATION
                | Chem.SanitizeFlags.SANITIZE_SETHYBRIDIZATION
                | Chem.SanitizeFlags.SANITIZE_SYMMRINGS,
                catchErrors=True,
            )
            mol.SetProp("partiallySanitized", "True")

        return mol

    @staticmethod
    def _get_mol(mol: Mol) -> Mol:
        """
        Get mol from SMILES (MCS seems to fail for equivalent unaligned molecules otherwise) unless maintaining
        partial sanitization required.
        """
        if "partiallySanitized" not in list(mol.GetPropNames()):
            return Chem.MolFromSmiles(Chem.MolToSmiles(mol))
        else:
            return mol

    def _calc_mcs(self, fragment: Mol, merge: Mol) -> Mol:
        """Function checks if there is a MCS between fragment A, B and the merge.
        Returns the molecule if there is and returns None if not.

        :param fragment:
        :param merge:

        :return: the MCS molecule or None
        :rtype: RDKit molecule or None
        """
        mcs = rdFMCS.FindMCS(
            [self._get_mol(fragment), self._get_mol(merge)], completeRingsOnly=True
        )
        mcs_mol = Chem.MolFromSmarts(mcs.smartsString)
        mcs_mol = self._sanitization(mcs_mol)
        return mcs_mol

    @staticmethod
    def _atom_remover(mol, pattern):
        matches = mol.GetSubstructMatches(pattern)
        if not matches:
            yield Chem.Mol(mol)
        for match in matches:
            res = Chem.RWMol(mol)
            res.BeginBatchEdit()
            for aid in match:
                res.RemoveAtom(aid)
            res.CommitBatchEdit()
            yield res

    @staticmethod
    def _check_synthon_mcs(
        mcs: Mol, min_atoms: int = config_filter.N_MCS_ATOMS
    ) -> bool:
        """
        Check if any non-carbon atoms OR if atom count >3.

        :param mcs: mcs mol
        :type mcs: RDKit molecule

        :return: whether molecule passes
        :rtype: bool
        """
        num_atoms = Lipinski.HeavyAtomCount(mcs)
        if num_atoms >= min_atoms:
            return True
        else:
            return False

    def filter_smi(
        self,
        smiles: str,
        synthon: str,
        fragmentA: Mol,
        fragmentB: Mol,
    ) -> bool:
        """
        Checks if molecule resemble expansions of fragment A (where fragment B is not contributing
        anything unique to the merge). Calculates the MCS between fragment A and the merge. Looks at the remainder
        of the merge and checks whether the MCS with the synthon is sensible (so fragment B has made a contribution -
        here if atom count is at least 3).

        :param smiles: smiles string of the merge
        :type smiles: str
        :param synthon: smiles string of the synthon
        :type synthon: str
        :param fragmentA: fragment A molecule
        :type fragmentA: RDKit molecule
        :param fragmentB: fragment B molecule
        :type fragmentB: RDKit molecule

        :return: whether molecule passes (True) or fails (False) filter
        :rtype: bool
        """
        merge = Chem.MolFromSmiles(smiles)
        mcs = self._calc_mcs(fragmentA, merge)
        mcs_removed_mols = [x for x in self._atom_remover(merge, mcs)]
        synthon = remove_xe(Chem.MolFromSmiles(synthon))
        result = False
        for _mol in mcs_removed_mols:
            mol = self._sanitization(_mol)
            synthon_mcs = self._calc_mcs(mol, synthon)
            if synthon_mcs:
                result = self._check_synthon_mcs(synthon_mcs)

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
            delayed(self.filter_smi)(smi, synthon, self._fragmentA, self._fragmentB)
            for smi, synthon in zip(self.smis, self.synthons)
        )
        return self.results, self.mols
