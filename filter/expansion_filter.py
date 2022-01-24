"""
Used to filter out compounds that look like elaborations rather than merges.
"""

from joblib import Parallel, delayed
from rdkit import Chem
from rdkit.Chem import AllChem, rdFMCS, rdShapeHelpers
from typing import Tuple

from filter.config_filter import config_filter
from filter.embedding_filter import add_coordinates, remove_xe
from filter.generic_filter import Filter_generic

# This function is currently a bit of a mess - there is a lot of exception handling involved that
# hasn't been sorted out yet


class ExpansionFilter(Filter_generic):

    def __init__(self, smis: list, synthons: list, fragmentA, fragmentB, proteinA=None, proteinB=None,
                 merge=None, mols=None, names=None):
        super().__init__(smis, synthons, fragmentA, fragmentB, proteinA, proteinB, merge, mols, names)
        self.results = None

    def _check_for_mcs(self, _mols: list):
        """Function checks if there is a MCS between fragment A, B and the merge.
        Returns the molecule if there is and returns None if not.

        :param _mols: list of mols (fragments and merge)
        :type _mols: list of RDKit molecules

        :return: the MCS molecule or None
        :rtype: RDKit molecule or None
        """
        mcs = rdFMCS.FindMCS(_mols, completeRingsOnly=True)
        mcs_mol = Chem.MolFromSmarts(mcs.smartsString)
        smi = Chem.MolToSmiles(mcs_mol)
        if smi == '':
            return None
        else:
            mcs_mol = Chem.MolFromSmiles(Chem.MolToSmiles(mcs_mol))
            return mcs_mol

    def filter_smi(self, merge: str, fragmentA, fragmentB, synthon: str, vol: float = config_filter.FRAGA_VOLUME)\
            -> bool:
        """
        Function ensures that merges that are just expansions of one fragment are removed.
        Removes part of the merge common to both fragments (by using MCS and checking overlap).
        Then calculates the proportion of the volume of the merge that comes from expanded fragment A.
        If more than 90%, these merges are ruled out.

        :param merge: merge SMILES string
        :type merge: str
        :param fragmentA: fragment A molecule
        :type fragmentA: RDKit molecule
        :param fragmentB: fragment B molecule
        :type fragmentB: RDKit molecule
        :param synthon: synthon SMILES string
        :type synthon: str

        :return: whether molecule passes (True) or fails (False) filter
        :rtype: bool
        """
        merge = Chem.MolFromSmiles(merge)
        synthon = Chem.MolFromSmiles(synthon)
        synthon = remove_xe(synthon)  # remove Xe from synthon

        try:
            # in few cases there is > 1 synthon match, just pass the molecule
            if len(merge.GetSubstructMatches(synthon)) > 1:
                result = True

            else:
                # check if there is a MCS (if not returns none)
                mols = [merge, fragmentA, fragmentB]
                mcs_mol = self._check_for_mcs(mols)

                if not mcs_mol:
                    # if no MCS, calculate volume and ratio
                    AllChem.EmbedMolecule(merge)
                    merge_vol = AllChem.ComputeMolVolume(merge)
                    fA_part = AllChem.DeleteSubstructs(merge, synthon)
                    fA_part_vol = AllChem.ComputeMolVolume(fA_part)

                    ratio = fA_part_vol / merge_vol

                    if ratio >= vol:
                        result = False
                    else:
                        result = True

                else:
                    AllChem.EmbedMolecule(merge)  # generate conformation of the merge

                    # check if there is overlap of the MCS in fragment A and B
                    mcs_A = add_coordinates(fragmentA, mcs_mol)
                    mcs_B = add_coordinates(fragmentB, mcs_mol)

                    dist = rdShapeHelpers.ShapeProtrudeDist(mcs_A, mcs_B)

                    if dist <= 0.5:  # if there is overlap, delete the substructure
                        merge_no_mcs = Chem.RWMol(merge)
                        mcs_matches = merge_no_mcs.GetSubstructMatches(mcs_mol)
                        for index in sorted(mcs_matches[0], reverse=True):
                            merge_no_mcs.RemoveAtom(index)

                        merge_vol = AllChem.ComputeMolVolume(merge_no_mcs)

                        # check if synthon is still there
                        if len(merge_no_mcs.GetSubstructMatches(synthon)) == 1:
                            # delete synthon and calculate the ratio
                            fA_part = AllChem.DeleteSubstructs(merge_no_mcs, synthon)
                            fA_part_vol = AllChem.ComputeMolVolume(fA_part)
                            ratio = fA_part_vol / merge_vol

                            if ratio >= vol:
                                result = False
                            else:
                                result = True

                        elif len(merge_no_mcs.GetSubstructMatches(synthon)) > 1:
                            # if >1 synthon substructure present, remove just one
                            fA_part = Chem.RWMol(merge_no_mcs)  # need editable mol
                            matches = fA_part.GetSubstructMatches(synthon)
                            for index in sorted(matches[0], reverse=True):
                                fA_part.RemoveAtom(index)

                            fA_part_vol = AllChem.ComputeMolVolume(fA_part)
                            ratio = fA_part_vol / merge_vol

                            if ratio >= vol:
                                result = False
                            else:
                                result = True

                        else:
                            result = False

                    else:   # if no overlap, calculate ratio of volumes
                        merge_vol = AllChem.ComputeMolVolume(merge)
                        fA_part = AllChem.DeleteSubstructs(merge, synthon)
                        fA_part_vol = AllChem.ComputeMolVolume(fA_part)

                        ratio = fA_part_vol / merge_vol

                        if ratio >= vol:
                            result = False
                        else:
                            result = True
        except:
            result = False

        return result

    def filter_all(self, cpus: int = config_filter.N_CPUS_FILTER_PAIR) -> Tuple[list, None]:
        """
        Runs the expansion filter on all the SMILES in parallel.

        :param cpus: number of CPUs for parallelization
        :type cpus: int

        :return: list of results (True or False); list of mols (None)
        :rtype: tuple
        """
        print('debugging', self._fragmentA)
        self.results = Parallel(n_jobs=cpus, backend='multiprocessing') \
            (delayed(self.filter_smi)(smi, self._fragmentA, self._fragmentB, synthon) for smi, synthon in
             zip(self.smis, self.synthons))
        return self.results, self.mols
