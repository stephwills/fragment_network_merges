"""
Used for filtering merges according to calculated descriptors.
"""

import argparse
import time

from dm_job_utilities.dm_log import DmLog
from filter.config_filter import config_filter
from filter.generic_scoring import Score_generic
from joblib import Parallel, delayed
from numpy import array
from oddt import fingerprints, toolkit
from oddt.toolkits.ob import Molecule
from rdkit import Chem


class IfpScore(Score_generic):
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
        mol_files=None,
        holo_files=None,
        apo_files=None,
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
            mol_files,
            holo_files,
            apo_files,
        )
        self.scores = None

    @staticmethod
    def _get_protein(protein: str) -> Molecule:
        """
        Function loads the protein from the pdb file.
        """
        protein = next(toolkit.readfile("pdb", protein))
        protein.protein = True
        return protein

    @staticmethod
    def _get_mol(mol: str) -> Molecule:
        """
        Function loads the molecule from the mol file.
        """
        return next(toolkit.readfile("mol", mol))

    @staticmethod
    def make_fp(protein: Molecule, mol: Molecule) -> array:
        """
        Function creates interaction fingerprint between the molecule and protein.

        :param protein: protein loaded from file
        :type protein: oddt.toolkits.mb.Molecule
        :param mol: mol loaded from file
        :type mol: oddt.toolkits.mb.Molecule

        :return: fingerprint
        :rtype: numpy.array
        """
        fp = fingerprints.InteractionFingerprint(mol, protein)
        return fp

    @staticmethod
    def calc_bond_percentage(fragment_fp: array, merge_fp: array) -> float:
        """
        Alternative to Tversky (which requires the fp to be binary).
        Calculate the percentage of the bonds made by the fragment that are preserved by
        the merge. Only checks bits in which the fragment makes an interaction but still
        accounts for the fact that this is a count vector.

        :param fragment_fp: interaction fingerprint of fragment
        :type fp: numpy.array
        :param merge_fp: interaction fingerprint of merge
        :type merge_fp: numpy.array

        :return: proportion of interactions maintained by the merge
        :rtype: float
        """
        frag_count = 0
        merge_count = 0

        for i, j in zip(fragment_fp, merge_fp):
            if i > 0:  # only check bits in which the fragment makes an interaction
                frag_count += i
                if i >= j:
                    merge_count += j
                else:
                    merge_count += i  # avoid percentage >100%

        if frag_count != 0:
            perc = merge_count / frag_count
        else:
            perc = 1

        return perc

    def ifp_score_avg(
        self, merge: str, fragmentA: str, fragmentB: str, prot: str
    ) -> float:
        """
        Function to calculate the average percentage bonds preserved between the merge and the fragment
        fingerprints. Calculates for each fragment and returns the average.

        :param merge: the file of the merge to filter
        :type merge: str
        :param fragmentA: the fragment A file
        :type fragmentA: str
        :param fragmentB: the fragment B file
        :type fragmentB: str
        :param prot: protein file
        :type prot: str

        :return: average percentage preserved
        :rtype: float
        """
        # load the molecules
        protein = self._get_protein(prot)
        merge_mol = self._get_mol(merge)
        fA_mol = self._get_mol(fragmentA)
        fB_mol = self._get_mol(fragmentB)
        # create all the fingerprints
        merge_fp, fA_fp, fB_fp = (
            self.make_fp(protein, merge_mol),
            self.make_fp(protein, fA_mol),
            self.make_fp(protein, fB_mol),
        )

        perc_A = self.calc_bond_percentage(fA_fp, merge_fp)
        perc_B = self.calc_bond_percentage(fB_fp, merge_fp)
        mean = (perc_A + perc_B) / 2

        return mean

    def score_mol(
        self,
        merge: str,
        fragmentA: str,
        fragmentB: str,
        prot: str,
    ) -> float:
        """
        Calculates the percentage of interactions made by the original fragments that are
        maintained by the merge. Can either calculate as the average of the score for each
        fragment or as the percentage of all interactions made.

        :param merge: the file of the merge to filter
        :type merge: str
        :param fragmentA: the fragment A file
        :type fragmentA: str
        :param fragmentB: the fragment B file
        :type fragmentB: str
        :param prot: protein file
        :type prot: str

        :return: percentage total bonds preserved
        :rtype: float
        """
        score = self.ifp_score_avg(merge, fragmentA, fragmentB, prot)

        return score

    def score_all(self, cpus: int = config_filter.N_CPUS_FILTER_PAIR) -> list:
        if not self.apo_files:
            # if no apo files, remove ligand from holo files to create apo files
            print("Removing ligands")
            self.get_apo_files()

        self.scores = Parallel(n_jobs=cpus, backend="multiprocessing")(
            delayed(self.score_mol)(mol_file, self.fragmentA, self.fragmentB, apo_file)
            for mol_file, apo_file in zip(self.mol_files, self.apo_files)
        )

        return self.scores


def main():
    parser = argparse.ArgumentParser(
        epilog="""
    python filter/ifp_score.py --input_file data/toFilter.sdf --output_file results.sdf
    --fragmentA_file nsp13-x0176_0B.mol --fragmentB_file nsp13-x0034_0B.mol --score_threshold 0.5
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
        "-a", "--fragmentA_file", required=True, help="fragment A mol file"
    )
    parser.add_argument(
        "-b", "--fragmentB_file", required=True, help="fragment B mol file"
    )
    parser.add_argument(
        "-t",
        "--score_threshold",
        type=float,
        default=None,
        help="minimum IFP score required for molecules to pass filter",
    )

    args = parser.parse_args()

    scorer = IfpScore()
    DmLog.emit_event("ifp_filter: ", args)

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
                    try:
                        minimised_mol_file = mol.GetProp("minimised_mol_file")
                        apo_file = mol.GetProp("apo_file")
                        score = scorer.score_mol(
                            minimised_mol_file,
                            args.fragmentA_file,
                            args.fragmentB_file,
                            apo_file,
                        )
                        if args.score_threshold:
                            if score >= args.score_threshold:
                                hits += 1
                                mol.SetProp("IFP score", score)
                                w.write(mol)
                        else:
                            hits += 1
                            mol.SetProp("IFP score", score)
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
