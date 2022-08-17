"""
Used for filtering merges according to calculated descriptors.
"""

import json
import os
import shutil

import numpy as np
from filter.config_filter import config_filter
from filter.generic_scoring import Score_generic
from joblib import Parallel, delayed
from rdkit import Chem
from scipy import stats
from utils.utils import get_distance


REACTION_INFO = {
    "Amidation": {
        "rxn": "[#6:1](=[#8:2])-[#8].[#7;H3,H2,H1:3]>>[#6:1](=[#8:2])-[#7:3]",
        "reactant1": Chem.MolFromSmiles("C(=O)O"),
        "reactant2": Chem.MolFromSmiles("CN"),
        "recursive_smarts1": Chem.MolFromSmarts("[O&$(OC(=O))]"),
        "recursive_smarts2": Chem.MolFromSmarts("[#7;H3,H2,H1:3]"),
    },
    "Amide_schotten-baumann": {
        "rxn": "[#7;H2,H1:3].[#6:1](=[#8:2])-[#17]>>[#6:1](=[#8:2])-[#7:3]",
        "reactant1": Chem.MolFromSmiles("CN"),
        "reactant2": Chem.MolFromSmiles("C(=O)Cl"),
        "recursive_smarts1": Chem.MolFromSmarts("[#7;H3,H2,H1:3]"),
        "recursive_smarts2": Chem.MolFromSmarts("[Cl&$(ClC(=O))]"),
    },
    "Reductive_amination": {
        "rxn": "[#6:2](=[#8])(-[#6:1]).[#7;H3,H2,H1:3]>>[#6:2](-[#6:1])-[#7:3]",
        "reactant1": Chem.MolFromSmiles("CC(=O)O"),
        "reactant2": Chem.MolFromSmiles("CN"),
        "recursive_smarts1": Chem.MolFromSmarts("[O&$(O=C(C))]"),
        "recursive_smarts2": Chem.MolFromSmarts("[#7;H3,H2,H1:3]"),
    },
    "N-nucleophilic_aromatic_substitution": {
        "rxn": "[#6:3]-[#7;H3,H2,H1:2].[c:1]-[F,Cl,Br,I]>>[#6:3]-[#7:2]-[c:1]",
        "reactant1": Chem.MolFromSmiles("CN"),
        "reactant2": Chem.MolFromSmiles("Fc1ccccc1"),
        "recursive_smarts1": Chem.MolFromSmarts("[N&$([#7;H3,H2,H1:2]-[#6:3])]"),
        "recursive_smarts2": Chem.MolFromSmarts("[F,Br,I,Cl&$([F,Br,I,Cl]c)]"),
    },
    "Sp2-sp2_Suzuki_coupling": {
        "rxn": "[c:1]-[F,Cl,Br,I].[#6:2]-[B]>>[c:1]-[#6:2]",
        "reactant1": Chem.MolFromSmiles("Fc1ccccc1"),
        "reactant2": Chem.MolFromSmiles("CB"),
        "recursive_smarts1": Chem.MolFromSmarts("[F,Br,I,Cl&$([F,Br,I,Cl]c)]"),
        "recursive_smarts2": Chem.MolFromSmarts("[B&$(Bc)]"),
    },
}


class ElaboratabilityScore(Score_generic):
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
    def get_attachment_points(mol, reaction_dict=REACTION_INFO):
        """
        Returns the indices of attachment point atoms by checking for substructure match with the reactants for some
        specified reactions. The index of the exact attachment point atom is identified by using recursive SMARTS (these
        have been written out manually).
        """
        attach_idxs = set()
        patts = ["recursive_smarts1", "recursive_smarts2"]
        for reaction in reaction_dict:
            for patt in patts:
                matches = mol.GetSubstructMatches(reaction_dict[reaction][patt])
                if len(matches) > 0:
                    # print(reaction_dict[reaction][patt], matches)
                    for match in matches:
                        attach_idxs.add(match[0])
        return list(attach_idxs)

    @staticmethod
    def calc_angle_between_points(attach_coords, atom1_coords, atom2_coords):
        """
        Calculate the angle between three atoms (two elaboration vectors pointed at two different protein atoms) using
        their coordinates.
        """
        a = np.array(atom1_coords)
        b = np.array(attach_coords)
        c = np.array(atom2_coords)

        ba = a - b
        bc = c - b

        cos_angle = np.dot(ba, bc) / (np.linalg.norm(ba) * np.linalg.norm(bc))
        angle = np.arccos(cos_angle)
        return np.degrees(angle)

    def dist_to_closest_spread_atoms(
        self,
        mol,
        atom_idx,
        prot,
        n_dists=config_filter.N_ELAB_DISTS,
        min_angle=config_filter.MIN_ELAB_ANGLE,
    ):
        conf = mol.GetConformer()
        pos = np.array(
            conf.GetAtomPosition(atom_idx)
        )  # get coordinates of the ligand exit vector atom

        dists = []
        prot_coords = []

        prot_conf = prot.GetConformer()
        for i in range(prot_conf.GetNumAtoms()):
            prot_pos = np.array(
                prot_conf.GetAtomPosition(i)
            )  # get the positions of all protein atoms
            prot_coords.append(prot_pos)
            dist = get_distance(
                pos, prot_pos
            )  # calculate the distance away from the ligand atom
            dists.append(dist)

        # take the geometric mean of the closest distances
        prot_coords = [
            coord for _, coord in sorted(zip(dists, prot_coords), key=lambda x: x[0])
        ]
        dists.sort()

        # record the closest 'spread out' dists (and the coords of the protein atoms)
        spread_dists = []
        spread_coords = []

        # keep iterating over the closest dists to get those that are not in the same direction
        c = 0

        for coord, dist in zip(prot_coords, dists):
            if len(spread_dists) == 0:
                spread_dists.append(dist)
                spread_coords.append(coord)
                c += 1

            else:
                if c < n_dists:
                    angles = []
                    for spread_dist, spread_coord in zip(spread_dists, spread_coords):
                        new_coord = coord
                        angle = self.calc_angle_between_points(
                            pos, spread_coord, new_coord
                        )
                        angles.append(angle)

                    # if the min angle against the already recorded atoms/dists is < the min angle, then add to the list
                    # print(angles, dist)
                    if min(angles) >= min_angle:
                        spread_dists.append(dist)
                        spread_coords.append(coord)
                        c += 1

        mean = stats.mstats.gmean(spread_dists)
        return mean

    def calculate_dists_to_nearest_atoms(
        self,
        mol,
        apo_file,
        attach_idxs,
        n_dists=config_filter.N_ELAB_DISTS,
        min_angle=config_filter.MIN_ELAB_ANGLE,
    ):
        prot = Chem.MolFromPDBFile(apo_file)
        dists = []
        for idx in attach_idxs:
            dists.append(
                self.dist_to_closest_spread_atoms(mol, idx, prot, n_dists, min_angle)
            )

        return dists

    def _write_and_copy_file(self, name, attach_idxs, dists, mol_file):
        d = {"attach_idxs": attach_idxs, "dists": dists}
        json_file = mol_file.replace(".minimised.mol", "_elab.json")
        if not os.path.exists(json_file):
            with open(json_file, "w") as f:
                json.dump(d, f)

        output_subdir = os.path.join(self.out_pair_dir, name.replace('_', '-'))
        if not os.path.exists(output_subdir):
            os.mkdir(output_subdir)

        outputdir_file = os.path.join(output_subdir, os.path.basename(json_file))
        shutil.copy(json_file, outputdir_file)

    def score_mol(
        self,
        name,
        mol,
        mol_file,
        apo_file,
        n_dists=config_filter.N_ELAB_DISTS,
        min_angle=config_filter.MIN_ELAB_ANGLE,
    ) -> float:
        """ """
        attach_idxs = self.get_attachment_points(mol)
        try:
            dists = self.calculate_dists_to_nearest_atoms(
                mol, apo_file, attach_idxs, n_dists, min_angle
            )
            self._write_and_copy_file(name, attach_idxs, dists, mol_file)
        except Exception as e:
            print('Failed for', name)
            print(e)
            print('Did not calc dists')
        return len(attach_idxs)

    def score_all(self, cpus: int = config_filter.N_CPUS_FILTER_PAIR) -> list:
        if not self.apo_files:
            # if no apo files, remove ligand from holo files to create apo files
            print("Removing ligands")
            self.get_apo_files()
            self.move_apo_files()

        self.scores = Parallel(n_jobs=cpus, backend="multiprocessing")(
            delayed(self.score_mol)(name, mol, mol_file, apo_file)
            for name, mol, mol_file, apo_file in zip(
                self.names, self.mols, self.mol_files, self.apo_files
            )
        )

        return self.scores
