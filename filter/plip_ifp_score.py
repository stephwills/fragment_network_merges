"""
Used for filtering merges according to calculated descriptors.
"""

import argparse
import json
import os
import shutil
import sys

from filter.config_filter import config_filter
from filter.generic_scoring import Score_generic
from joblib import Parallel, delayed
from plip.structure.preparation import PDBComplex
from pymol import cmd
from utils.utils import load_json

ALL_INTERACTIONS = [
    "hydrophobic_contacts",
    "pi_stacking",
    "pi_cation",
    "hbond_pdon",
    "hbond_ldon",
    "saltbridge_lneg",
    "saltbridge_pneg",
    "saltbridge_metal",
    "halogen",
]


class PlipIfpScore(Score_generic):
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
    def make_complex_file(fragment: str, fragment_file: str, apo_file: str) -> str:
        """
        Create pdb file containing complex of fragment and the apo minimized protein for PLIP
        """
        new_filename = apo_file.replace("nolig.pdb", f"f{fragment}.pdb")
        if not os.path.exists(new_filename):
            cmd.reinitialize()
            cmd.load(apo_file, "prot")
            cmd.load(fragment_file, "ligand")
            cmd.create("complex", "ligand, prot")
            cmd.save(new_filename, "complex")
        return new_filename

    @staticmethod
    def get_interaction_set(complex_file):
        """
        Load the PDBComplex and get the interaction set to process.
        """
        mol = PDBComplex()
        mol.load_pdb(complex_file)
        mol.analyze()
        bsid = [":".join([x.hetid, x.chain, str(x.position)]) for x in mol.ligands][
            0
        ]  # get binding site id
        interactions = mol.interaction_sets[bsid]
        return interactions

    def get_interactions(self, complex_file):
        """
        Load the complex pdb file using PLIP and analyse to get interactions in the form a dictionary, with interaction
        type as the key and values as list of residues that make that interaction
        """
        interactions = self.get_interaction_set(complex_file)

        hydrophobic_contacts = [
            f"{int_type.restype}-{int_type.resnr}"
            for int_type in interactions.hydrophobic_contacts
        ]
        pi_stacking = [
            f"{int_type.restype}-{int_type.resnr}"
            for int_type in interactions.pistacking
        ]
        pi_cation = [
            f"{int_type.restype}-{int_type.resnr}"
            for int_type in interactions.pication_laro
        ] + [
            f"{int_type.restype}-{int_type.resnr}"
            for int_type in interactions.pication_paro
        ]
        hbond_pdon = [
            f"{int_type.restype}-{int_type.resnr}"
            for int_type in interactions.hbonds_pdon
        ]
        hbond_ldon = [
            f"{int_type.restype}-{int_type.resnr}"
            for int_type in interactions.hbonds_ldon
        ]
        saltbridge_lneg = [
            f"{int_type.restype}-{int_type.resnr}"
            for int_type in interactions.saltbridge_lneg
        ]
        saltbridge_pneg = [
            f"{int_type.restype}-{int_type.resnr}"
            for int_type in interactions.saltbridge_pneg
        ]
        saltbridge_metal = [
            f"{int_type.restype}-{int_type.resnr}"
            for int_type in interactions.metal_complexes
        ]
        halogen = [
            f"{int_type.restype}-{int_type.resnr}"
            for int_type in interactions.halogen_bonds
        ]

        interaction_dict = {
            "hydrophobic_contacts": hydrophobic_contacts,
            "pi_stacking": pi_stacking,
            "pi_cation": pi_cation,
            "hbond_pdon": hbond_pdon,
            "hbond_ldon": hbond_ldon,
            "saltbridge_lneg": saltbridge_lneg,
            "saltbridge_pneg": saltbridge_pneg,
            "saltbridge_metal": saltbridge_metal,
            "halogen": halogen,
        }

        return interaction_dict

    def write_interaction_file(
        self, name, fragmentA_file, fragmentB_file, holo_file, apo_file
    ):
        json_file = holo_file.replace(".holo_minimised.pdb", "_interactions.json")
        if not os.path.exists(json_file):
            fragmentA_complex_file = self.make_complex_file(
                "A", fragmentA_file, apo_file
            )
            fragmentB_complex_file = self.make_complex_file(
                "B", fragmentB_file, apo_file
            )

            fragA_interactions = self.get_interactions(fragmentA_complex_file)
            fragB_interactions = self.get_interactions(fragmentB_complex_file)
            merge_interactions = self.get_interactions(holo_file)

            d = {
                "merge": merge_interactions,
                "fragA": fragA_interactions,
                "fragB": fragB_interactions,
            }

            with open(json_file, "w") as f:
                json.dump(d, f)

            self._copy_files(
                name, [json_file, fragmentA_complex_file, fragmentB_complex_file]
            )

        return json_file

    def _copy_files(self, name: str, workdir_files):
        """
        Move the interaction and complex files to the output_dir
        """
        output_subdir = os.path.join(self.out_pair_dir, name.replace("_", "-"))
        if not os.path.exists(output_subdir):
            os.mkdir(output_subdir)

        outputdir_files = [
            os.path.join(output_subdir, os.path.basename(f)) for f in workdir_files
        ]
        for workdir_file, outputdir_file in zip(workdir_files, outputdir_files):
            shutil.copy(workdir_file, outputdir_file)

    def ifp_scorer(self, interactions_file):
        data = load_json(interactions_file)

        proportions = []

        for fragment in ["fragA", "fragB"]:
            fragment_total = 0
            merge_total = 0
            for int_type in ALL_INTERACTIONS:

                residues = data[fragment][int_type]
                merge_residues = data["merge"][int_type]

                if len(residues) > 0:
                    for residue in set(residues):
                        frag_count = residues.count(residue)
                        merge_count = merge_residues.count(residue)

                        if frag_count < merge_count:
                            fragment_total += frag_count
                            merge_total += frag_count
                        else:
                            fragment_total += frag_count
                            merge_total += merge_count

            if fragment_total == 0:
                # if the fragment doesn't make any interactions, set proportion as 0.5
                # not a perfect fix - TODO: change this
                proportion = 0.5
            else:
                proportion = merge_total / fragment_total
            proportions.append(proportion)

        mean = (proportions[0] * proportions[1]) ** 0.5
        return mean

    def score_mol(
        self,
        name: str,
        holo_file: str,
        apo_file: str,
        fragmentA: str,
        fragmentB: str,
    ) -> float:
        """ """
        interaction_file = self.write_interaction_file(
            name, fragmentA, fragmentB, holo_file, apo_file
        )
        score = self.ifp_scorer(interaction_file)

        return score

    def score_all(self, cpus: int = config_filter.N_CPUS_FILTER_PAIR) -> list:
        if not self.apo_files:
            # if no apo files, remove ligand from holo files to create apo files
            print("Removing ligands")
            self.get_apo_files()
            self.move_apo_files()

        self.scores = Parallel(n_jobs=cpus, backend="multiprocessing")(
            delayed(self.score_mol)(
                name, holo_file, apo_file, self.fragmentA, self.fragmentB
            )
            for name, holo_file, apo_file in zip(
                self.names, self.holo_files, self.apo_files
            )
        )

        return self.scores

    def filter_mol(
        self,
        name: str,
        holo_file: str,
        apo_file: str,
        fragmentA: str,
        fragmentB: str,
        score_threshold: float,
    ) -> float:
        """
        Filter using the scoring function.
        """
        interaction_file = self.write_interaction_file(
            name, fragmentA, fragmentB, holo_file, apo_file
        )
        score = self.ifp_scorer(interaction_file)
        if score >= score_threshold:
            result = True
        else:
            result = False
        return result

    def filter_all(self, cpus: int = config_filter.N_CPUS_FILTER_PAIR, *args):
        """
        Runs the PLIP IFP score filter on all molecules in parallel.
        """
        if not self.apo_files:
            # if no apo files, remove ligand from holo files to create apo files
            print("Removing ligands")
            self.get_apo_files()
            # self.move_apo_files()

        results = Parallel(n_jobs=cpus, backend="multiprocessing")(
            delayed(self.filter_mol)(
                name, holo_file, apo_file, self.fragmentA, self.fragmentB, *args
            )
            for name, holo_file, apo_file in zip(
                self.names, self.holo_files, self.apo_files
            )
        )

        return results, self.mols


def parse_args(args):
    """
    Parse command line arguments
    """
    parser = argparse.ArgumentParser(
        epilog="""
    python filter/plip_ifp_score.py --input_file data/toFilter.sdf --output_file results.sdf --score_threshold 0.6
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
        "-t",
        "--score_threshold",
        required=False,
        help="mean proportion of interactions that need to be maintained to pass the filter",
        type=float,
        default=0.5,
    )
    parser.add_argument(
        "-W", "--working_dir", required=True, help="where to save intermediate files"
    )
    parser.add_argument(
        "-O", "--output_dir", required=True, help="where to save the final files"
    )
    return parser.parse_args(args)


def main():
    from filter.generic_squonk import Squonk_generic

    args = parse_args(sys.argv[1:])
    job = Squonk_generic(
        "PlipIfpScore",
        args,
        args.input_file,
        args.output_file,
        args.fragmentA_file,
        args.fragmentB_file,
        work_pair_dir=args.working_dir,
        out_pair_dir=args.output_dir,
        score=True,
    )
    job.execute_job(args.n_cpus, args.score_threshold)


if __name__ == "__main__":
    main()
