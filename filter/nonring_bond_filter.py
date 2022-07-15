"""
For filtering molecules with many consecutive non-ring bonds atoms.
"""

import argparse
import sys
import time
from typing import Tuple

from dm_job_utilities.dm_log import DmLog
from filter.config_filter import config_filter
from filter.generic_filter import Filter_generic
from joblib import Parallel, delayed
from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem, rdmolops, Mol
from rdkit.Chem.rdMolDescriptors import CalcNumHeavyAtoms, CalcNumRings

# uses code from
# https://iwatobipen.wordpress.com/2020/01/23/cut-molecule-to-ring-and-linker-with-rdkit-rdkit-chemoinformatics-memo/
RDLogger.DisableLog("rdApp.*")


class NonringBondFilter(Filter_generic):
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
    def is_in_same_ring(idx1: int, idx2: int, bond_rings) -> bool:
        """
        Checks if two atoms (using atom indices) are in the same ring
        """
        for bond_ring in bond_rings:
            if idx1 in bond_ring and idx2 in bond_ring:
                return True
        return False

    @staticmethod
    def get_fragments(mol: Mol, bonds: list) -> Tuple[list, list]:
        """
        Cleaves the molecule into fragments at specified bonds.
        """
        _fragments = Chem.FragmentOnBonds(mol, bonds)
        fragments = Chem.GetMolFrags(_fragments, asMols=True)
        linkers = []
        sidechains = []

        for fragment in fragments:
            # check the fragment does not contain rings
            numRings = CalcNumRings(fragment)
            if numRings == 0:
                smiles = Chem.MolToSmiles(fragment)
                # check number of attachment points (one for sidechains)
                if smiles.count("*") == 1:
                    sidechains.append(fragment)
                else:
                    linkers.append(fragment)

        return linkers, sidechains

    @staticmethod
    def max_path(fragment: Mol) -> int:
        """
        Get the maximum path length within the substructure.
        """
        maxPath = 0
        numAtoms = CalcNumHeavyAtoms(fragment)
        for i in range(numAtoms):  # look for paths up the length of number of atoms
            # get max consecutive path length in linker/sidechain
            paths = rdmolops.FindAllPathsOfLengthN(fragment, i, useBonds=True)
            if len(paths) > 0:
                maxPath = i

        return maxPath

    def get_bonds_to_cleave(self, mol: Mol):
        """
        Get the bonds that are attached to rings (cleaving at these bonds will give the linker and sidechain
        substructures.
        """
        cleave_bonds = []

        # record original indices of atoms and bonds as properties
        for atom in mol.GetAtoms():
            atom.SetIntProp("orig_idx", atom.GetIdx())
        for bond in mol.GetBonds():
            bond.SetIntProp("orig_idx", bond.GetIdx())

        # get set of bonds that are in rings
        ring_info = mol.GetRingInfo()
        bond_rings = (
            ring_info.BondRings()
        )  # gives tuples with bonds that are in the rings
        ring_bonds = set()
        for ring_bond_idxs in bond_rings:  # get set
            for idx in ring_bond_idxs:
                ring_bonds.add(idx)

        # get the non-ring bonds
        all_bonds_idx = [bond.GetIdx() for bond in mol.GetBonds()]
        non_ring_bonds = set(all_bonds_idx) - ring_bonds

        for bond_idx in non_ring_bonds:
            # get the indices of the atoms on each end of the bond
            bgn_idx = mol.GetBondWithIdx(bond_idx).GetBeginAtomIdx()
            end_idx = mol.GetBondWithIdx(bond_idx).GetEndAtomIdx()
            if (
                mol.GetBondWithIdx(bond_idx).GetBondTypeAsDouble() == 1.0
            ):  # if a single bond
                # GetBondTypeAsDouble() - 1.0 for single, 1.5 for aromatic, 2.0 for double
                # if one of the atoms on the side of the bond is in the ring, append to cleave_bonds
                if (
                    mol.GetAtomWithIdx(bgn_idx).IsInRing()
                    + mol.GetAtomWithIdx(end_idx).IsInRing()
                    == 1
                ):
                    bond = mol.GetBondWithIdx(bond_idx)
                    orig_idx = bond.GetIntProp("orig_idx")
                    cleave_bonds.append(orig_idx)
                # if the bond atoms are in two different rings, append to cleave_bonds
                elif (
                    not self.is_in_same_ring(bgn_idx, end_idx, bond_rings)
                    and mol.GetAtomWithIdx(bgn_idx).IsInRing()
                    + mol.GetAtomWithIdx(end_idx).IsInRing()
                    == 2
                ):
                    bond = mol.GetBondWithIdx(bond_idx)
                    orig_idx = bond.GetIntProp("orig_idx")
                    cleave_bonds.append(orig_idx)
        return cleave_bonds

    def filter_smi(
        self,
        smi: str,
        linker_path_threshold: int = config_filter.LINKER_PATH_THRESHOLD,
        sidechain_path_threshold: int = config_filter.SIDECHAIN_PATH_THRESHOLD,
    ):
        """
        Checks if the SMILES fails for number of consecutive non-ring atoms (different thresholds set dependent on
        whether considering linkers or sidechains.

        :param smi: smiles of the merge
        :type smi: str
        :param linker_path_threshold: max length of linkers (joined to two rings)
        :type linker_path_threshold: int
        :param sidechain_path_threshold: max length of sidechains (joined to one ring)
        :type sidechain_path_threshold: int

        :return: whether molecule passes (True) or fails (False) filter
        :rtype: bool
        """
        mol = Chem.MolFromSmiles(smi)
        try:
            # get the bonds to cleave for the molecule (joined to a ring)
            bonds = self.get_bonds_to_cleave(mol)
            if len(bonds) == 0:
                return True
            else:
                linkers, sidechains = self.get_fragments(mol, bonds)

                max_linker_path = 0
                max_sidechain_path = 0

                if len(linkers) > 0:
                    linker_paths = [self.max_path(linker) for linker in linkers]
                    if len(linker_paths) > 0:
                        max_linker_path = max(linker_paths)

                if len(sidechains) > 0:
                    sidechain_paths = [
                        self.max_path(sidechain) for sidechain in sidechains
                    ]
                    if len(sidechain_paths) > 0:
                        max_sidechain_path = max(sidechain_paths)

                if (
                    max_linker_path < linker_path_threshold
                    and max_sidechain_path < sidechain_path_threshold
                ):
                    return True

                else:
                    return False

        except:
            return True

    def filter_all(
        self, cpus: int = config_filter.N_CPUS_FILTER_PAIR
    ) -> Tuple[list, None]:
        """
        Runs the nonring_bond_filter on all the SMILES in parallel.

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
    python filter/nonring_bond_filter.py --input_file data/toFilter.sdf --output_file results.sdf
    --linker_path_threshold 8 --sidechain_path_threshold 6
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
        "-l",
        "--linker_path_threshold",
        default=8,
        help="threshold for length of linker substructures",
    )
    parser.add_argument(
        "-s",
        "--sidechain_path_threshold",
        default=6,
        help="threshold for length of sidechain substructures",
    )
    return parser.parse_args(args)


def main():
    args = parse_args(sys.argv[1:])
    filter = NonringBondFilter()
    DmLog.emit_event("nonring_bond_filter: ", args)

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
                    res = filter.filter_smi(
                        smi, args.linker_path_threshold, args.sidechain_path_threshold
                    )
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
