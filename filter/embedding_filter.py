"""
Used for 3D filtering of fragment merges by constrained embedding.
"""

import os
import sys
import argparse
import time
from typing import Tuple

import numpy as np
from filter.config_filter import config_filter
from filter.generic_filter import Filter_generic
from joblib import Parallel, delayed
from rdkit import Chem, RDLogger
from rdkit.Chem import Mol, rdForceFieldHelpers, rdmolfiles
from utils.filter_utils import (
    ConstrainedEmbedMatches,
    add_coordinates,
    calc_energy,
    calc_unconstrained_energy,
    get_mcs,
    remove_xe,
)
from utils.utils import get_distance

RDLogger.DisableLog("rdApp.*")


class EmbeddingFilter(Filter_generic):
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

    @staticmethod
    def check_overlap(
        molA: Mol, molB: Mol, atom_clash_dist: float = config_filter.ATOM_CLASH_DIST
    ) -> Tuple[Mol, Mol]:
        """
        Function checks if parts of two molecules overlap. If atoms overlap, then they are removed from one of the
        molecules.
        """
        A = Chem.RWMol(molA)
        B = Chem.RWMol(molB)
        confA = A.GetConformer()
        confB = B.GetConformer()
        clashes = []
        for i in range(confA.GetNumAtoms()):
            posA = np.array(confA.GetAtomPosition(i))
            for j in range(confB.GetNumAtoms()):
                posB = np.array(confB.GetAtomPosition(j))
                dist = get_distance(posA, posB)
                # if atoms closer than atom_clash_dist, recognised as clashes
                if dist < atom_clash_dist:
                    clashes.append(i)
                    break
        if clashes:  # remove clashing atoms from one of the molecules
            # sort in descending order (so atoms removed correctly)
            sorted_list = sorted(clashes, reverse=True)
            for clash in sorted_list:
                A.RemoveAtom(clash)  # remove atoms from one fragment
        return A, B

    def embedding(
        self,
        fragA: Mol,
        fragB: Mol,
        merge_mol: Mol,
        synth,
        atom_clash_dist: float = config_filter.ATOM_CLASH_DIST,
    ) -> list:
        """
        Function to embed the full molecule, constraining the atoms that came from each fragment.
        The atoms that came from each fragment are retrieved, and the 3D coordinates
        are added from the original structures.

        :param fragA: original fragment A with 3D conformation
        :type fragA: rdkit.Chem.Mol
        :param fragB: original fragment B with 3D conformation
        :type fragB: rdkit.Chem.Mol
        :param merge_mol: proposed merge
        :type merge_mol: rdkit.Chem.Mol
        :param synth: synthon from fragment B used for expansion if frag net; None if sim search
        :type synth: rdkit.Chem.Mol or None
        :param atom_clash_dist: minimum distance between atoms to be considered as clashing
        :type atom_clash_dist: float

        :return: list of embedded molecules (if embedding was successful, otherwise empty)
        :rtype: list
        """
        # substructure match for fragment A
        # identify the atoms that came from fragment A
        mcsA = get_mcs(merge_mol, fragA)
        # get all possible matches with fragment A
        mcsA_matches = fragA.GetSubstructMatches(mcsA)
        # for each, add coordinates
        mcsA_mols = [add_coordinates(fragA, mcsA, matches) for matches in mcsA_matches]

        if synth:
            # substructure match for fragment B
            synthB = remove_xe(synth)  # remove the xenon from the synthon
            synthB_matches = fragB.GetSubstructMatches(synthB)
            synthB_mols = [
                add_coordinates(fragB, synthB, matches) for matches in synthB_matches
            ]

            # get all possible ref coordinates
            ref_mols = []  # store ref molecules with different coordinates
            for _mcsA_mol in mcsA_mols:
                for _synthB_mol in synthB_mols:
                    # check if atoms overlap
                    mcsA_mol, synthB_mol = self.check_overlap(
                        _mcsA_mol, _synthB_mol, atom_clash_dist
                    )
                    # combine substructures to get ref molecule
                    ref_mol = Chem.CombineMols(mcsA_mol, synthB_mol)
                    ref_mols.append(ref_mol)

        else:
            # substructure match for fragment B
            # identify the atoms that came from fragment B
            mcsB = get_mcs(merge_mol, fragB)
            # get all possible matches with fragment B
            mcsB_matches = fragB.GetSubstructMatches(mcsB)
            # for each, add coordinates
            mcsB_mols = [
                add_coordinates(fragB, mcsB, matches) for matches in mcsB_matches
            ]

            # get all possible ref coordinates
            ref_mols = []  # store ref molecules with different coordinates
            for _mcsA_mol in mcsA_mols:
                for _mcsB_mol in mcsB_mols:
                    # check if atoms overlap
                    mcsA_mol, mcsB_mol = self.check_overlap(
                        _mcsA_mol, _mcsB_mol, atom_clash_dist
                    )
                    ref_mol = Chem.CombineMols(mcsA_mol, mcsB_mol)
                    ref_mols.append(ref_mol)

        # Get substructure matches for merge and embed with all sets of coordinates
        all_merge_matches = [
            merge_mol.GetSubstructMatches(ref_mol) for ref_mol in ref_mols
        ]
        embedded_mols = []
        n_embeddings = 0
        for merge_matches, ref_mol in zip(all_merge_matches, ref_mols):
            for matches in merge_matches:
                n_embeddings += 1
                merge = Chem.RWMol(merge_mol)
                # merge = Chem.AddHs(merge)
                try:
                    embedded_mol = ConstrainedEmbedMatches(
                        merge, ref_mol, matches, randomseed=42
                    )
                    rdForceFieldHelpers.MMFFOptimizeMolecule(embedded_mol)
                    # embedded_mol = Chem.RemoveHs(embedded_mol)
                    embedded_mols.append(embedded_mol)
                except ValueError:
                    pass

        # print(f"{len(embedded_mols)}/{n_embeddings} embeddings were successful")

        return embedded_mols

    def filter_smi(
        self,
        merge: str,
        fragA: Mol,
        fragB: Mol,
        synthon,  # str if frag net; None if sim search
        energy_threshold: float = config_filter.ENERGY_THRESHOLD,
        n_conf: int = config_filter.N_CONFORMATIONS,
        atom_clash_dist: float = config_filter.ATOM_CLASH_DIST,
    ) -> Tuple[bool, None]:
        """
        Runs the filter by embedding (if possible) the molecule, calculating the energy of the constrained
        conformation, and comparing with the avg energy of the unconstrained conformations (specified by n_conf). If
        the energy ratio is greater than a specified threshold for the constrained molecule, the molecule is filtered
        out.

        :param merge: merge SMILES string
        :type merge: str
        :param fragA: fragment A molecule
        :type fragA: rdkit.Chem.Mol
        :param fragB: fragment B molecule
        :type fragB: rdkit.Chem.Mol
        :param synthon: synthon SMILES string
        :type synthon: str
        :param energy_threshold: threshold of constrained energy versus unconstrained energy ratio
        :type energy_threshold: float
        :param n_conf: number of unconstrained conformations to generate
        :type n_conf: int
        :param atom_clash_dist: minimum distance between atoms to be considered as clashing
        :type atom_clash_dist: float

        :return: whether molecule passes (True) or fails (False) filter
        :rtype: bool
        """
        try:
            merge = Chem.MolFromSmiles(merge)
            if synthon:
                synthon = Chem.MolFromSmiles(synthon)
            embedded_mols = self.embedding(
                fragA, fragB, merge, synthon, atom_clash_dist
            )

            result, embedded = False, None
            if len(embedded_mols) == 0:
                result = False
                embedded = None

            else:
                for embedded_mol in embedded_mols:
                    # energy of constrained conformation
                    const_energy = calc_energy(embedded_mol)
                    # energy of avg unconstrained conformation
                    unconst_energy = calc_unconstrained_energy(merge, n_conf)
                    # if the energy of the constrained conformation is less, then pass filter
                    if const_energy <= unconst_energy:
                        result = True
                        embedded = embedded_mol
                        break
                    else:
                        # if constrained energy > energy-threshold-fold greater, then fail filter
                        ratio = const_energy / unconst_energy
                        if ratio >= energy_threshold:
                            result = False
                            embedded = None
                        else:
                            result = True
                            embedded = embedded_mol
                            break

            return result, embedded

        except Exception as e:  # pass molecules that break the filter
            print("Failed for merge", merge)
            print(e)
            return False, None

    def filter_all(
        self, cpus: int = config_filter.N_CPUS_FILTER_PAIR, *args
    ) -> Tuple[list, list]:
        """
        Runs the overlap filter on all the SMILES in parallel.

        :param cpus: number of CPUs for parallelization
        :type cpus: int

        :return: list of results (True or False); list of mols (RDKit molecules)
        :rtype: tuple
        """
        if self.synthons:
            res = Parallel(n_jobs=cpus, backend="multiprocessing")(
                delayed(self.filter_smi)(
                    smi, self._fragmentA, self._fragmentB, synthon, *args
                )
                for smi, synthon in zip(self.smis, self.synthons)
            )
        else:
            res = Parallel(n_jobs=cpus, backend="multiprocessing")(
                delayed(self.filter_smi)(
                    smi, self._fragmentA, self._fragmentB, None, *args
                )
                for smi in self.smis
            )
        self.results = [r[0] for r in res]
        self.mols = [r[1] for r in res]
        return self.results, self.mols


def parse_args(args):
    parser = argparse.ArgumentParser(
        epilog="""
        python filter/embedding_filter.py --input_file data/toFilter.sdf --output_file results.sdf
        --fragmentA_file fragmentA.mol --fragmentB_file fragmaentB.mol --energy_threshold 10
        --n_conformations 50 --atom_clash_dist 1.0
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
        type=int,
        default=10,
        help="threshold for fold difference in energy of unconstrained versus constrained conformations",
    )
    parser.add_argument(
        "-n",
        "--n_conformations",
        type=int,
        default=50,
        help="number of unconstrained conformations to generate",
    )
    parser.add_argument(
        "-d",
        "--atom_clash_dist",
        type=float,
        default=1.0,
        help="max distance between fragment A and B atoms to remove from ref mol for embedding",
    )
    return parser.parse_args(args)


def main():
    from filter.generic_squonk import Squonk_generic

    args = parse_args(sys.argv[1:])
    job = Squonk_generic(
        "EmbeddingFilter",
        args,
        args.input_file,
        args.output_file,
        args.fragmentA_file,
        args.fragmentB_file,
    )
    job.execute_job(
        args.n_cpus, args.energy_threshold, args.n_conformations, args.atom_clash_dist
    )


if __name__ == "__main__":
    main()
