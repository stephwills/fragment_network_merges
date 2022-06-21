"""
Used to filter out compounds that look like elaborations rather than merges.
"""

import argparse
import time
from typing import Tuple

from dm_job_utilities.dm_log import DmLog
from filter.config_filter import config_filter
from filter.generic_filter import Filter_generic
from joblib import Parallel, delayed
from rdkit import Chem
from rdkit.Chem import Lipinski, Mol, rdmolfiles
from utils.filter_utils import get_mcs, remove_xe


class ExpansionFilter(Filter_generic):
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
    def _get_mol(mol: Mol) -> Mol:
        """
        Get mol from SMILES (MCS seems to fail for equivalent unaligned molecules otherwise) unless maintaining
        partial sanitization required.
        """
        if "partiallySanitized" not in list(mol.GetPropNames()):
            return Chem.MolFromSmiles(Chem.MolToSmiles(mol))
        else:
            return mol

    @staticmethod
    def _atom_remover(mol, pattern):
        """
        Removes atoms from molecule that match pattern.
        """
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
        min_atoms: int = config_filter.N_MCS_ATOMS,
    ) -> bool:
        """
        Checks if molecule resemble expansions of fragment A (where fragment B is not contributing anything unique to
        the merge). Calculates the MCS between fragment A and the merge. Looks at the remainder of the merge and checks
        whether the MCS with the synthon is sensible (so fragment B has made a contribution - here if atom count is at
        least 3).

        :param smiles: smiles string of the merge
        :type smiles: str
        :param synthon: smiles string of the synthon
        :type synthon: str
        :param fragmentA: fragment A molecule
        :type fragmentA: RDKit molecule
        :param fragmentB: fragment B molecule
        :type fragmentB: RDKit molecule
        :param min_atoms: min number of atoms contributed from synthon
        :type min_atoms: int

        :return: whether molecule passes (True) or fails (False) filter
        :rtype: bool
        """
        merge = Chem.MolFromSmiles(smiles)
        mcs = get_mcs(fragmentA, merge)
        mcs_removed_mols = [x for x in self._atom_remover(merge, mcs)]
        synthon = remove_xe(Chem.MolFromSmiles(synthon))
        result = False
        for mol in mcs_removed_mols:
            synthon_mcs = get_mcs(mol, synthon)
            if synthon_mcs:
                result = self._check_synthon_mcs(synthon_mcs)

        return result

    def filter_all(
        self, cpus: int = config_filter.N_CPUS_FILTER_PAIR, **kwargs
    ) -> Tuple[list, None]:
        """
        Runs the expansion filter on all the SMILES in parallel.

        :param cpus: number of CPUs for parallelization
        :type cpus: int

        :return: list of results (True or False); list of mols (None)
        :rtype: tuple
        """
        self.results = Parallel(n_jobs=cpus, backend="multiprocessing")(
            delayed(self.filter_smi)(smi, synthon, self._fragmentA, self._fragmentB, **kwargs)
            for smi, synthon in zip(self.smis, self.synthons)
        )
        return self.results, self.mols


def main():
    parser = argparse.ArgumentParser(
        epilog="""
    python filter/expansion_filter.py --input_file data/toFilter.sdf --output_file results.sdf
    --fragmentA_file nsp13-x0176_0B.mol --min_atoms 3
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
        "-A", "--fragmentA_file", required=True, help="fragment A mol file"
    )
    parser.add_argument(
        "-B", "--fragmentB_file", required=True, help="fragment B mol file"
    )
    parser.add_argument(
        "-a",
        "--min_atoms",
        type=int,
        default=3,
        help="minimum atom contribution from fragment B",
    )

    args = parser.parse_args()

    filter = ExpansionFilter()
    DmLog.emit_event("expansion_filter: ", args)

    start = time.time()
    count = 0
    hits = 0
    errors = 0

    fragmentA = rdmolfiles.MolFromMolFile(args.fragmentA_file)
    fragmentB = rdmolfiles.MolFromMolFile(args.fragmentB_file)

    with Chem.SDWriter(args.output_file) as w:
        with Chem.SDMolSupplier(args.input_file) as suppl:
            for mol in suppl:
                if mol is None:
                    continue
                else:
                    count += 1
                    try:
                        smi = Chem.MolToSmiles(mol)
                        synthon = mol.GetProp("synthon")
                        res = filter.filter_smi(
                            smi, synthon, fragmentA, fragmentB, args.min_atoms
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
