"""
Used for 3D filtering of fragment merges by constrained embedding.
"""

import numpy as np

from joblib import Parallel, delayed
from rdkit import Chem
from rdkit.Chem.AllChem import *
from rdkit import RDLogger
from rdkit.Chem import AllChem, rdFMCS, rdForceFieldHelpers, Mol
from typing import Tuple

from filter.config_filter import config_filter
from filter.generic_filter import Filter_generic

RDLogger.DisableLog("rdApp.*")


def ConstrainedEmbedMatches(
    mol,
    core,
    match,
    useTethers=True,
    coreConfId=-1,
    randomseed=2342,
    getForceField=UFFGetMoleculeForceField,
    **kwargs,
):
    """
    Adapted from original RDKit function rdkit.Chem.AllChem.ConstrainedEmbed to iterate through multiple substructure
    matches
    """
    coordMap = {}
    coreConf = core.GetConformer(coreConfId)
    for i, idxI in enumerate(match):
        corePtI = coreConf.GetAtomPosition(i)
        coordMap[idxI] = corePtI

    ci = EmbedMolecule(mol, coordMap=coordMap, randomSeed=randomseed, **kwargs)
    if ci < 0:
        raise ValueError("Could not embed molecule.")

    algMap = [(j, i) for i, j in enumerate(match)]

    if not useTethers:
        # clean up the conformation
        ff = getForceField(mol, confId=0)
        for i, idxI in enumerate(match):
            for j in range(i + 1, len(match)):
                idxJ = match[j]
                d = coordMap[idxI].Distance(coordMap[idxJ])
                ff.AddDistanceConstraint(idxI, idxJ, d, d, 100.0)
        ff.Initialize()
        n = 4
        more = ff.Minimize()
        while more and n:
            more = ff.Minimize()
            n -= 1
        # rotate the embedded conformation onto the core:
        rms = AlignMol(mol, core, atomMap=algMap)
    else:
        # rotate the embedded conformation onto the core:
        rms = AlignMol(mol, core, atomMap=algMap)
        ff = getForceField(mol, confId=0)
        conf = core.GetConformer()
        for i in range(core.GetNumAtoms()):
            p = conf.GetAtomPosition(i)
            pIdx = ff.AddExtraPoint(p.x, p.y, p.z, fixed=True) - 1
            ff.AddDistanceConstraint(pIdx, match[i], 0, 0, 100.0)
        ff.Initialize()
        n = 4
        more = ff.Minimize(energyTol=1e-4, forceTol=1e-3)
        while more and n:
            more = ff.Minimize(energyTol=1e-4, forceTol=1e-3)
            n -= 1
        # realign
        rms = AlignMol(mol, core, atomMap=algMap)
    mol.SetProp("EmbedRMS", str(rms))
    return mol


def add_coordinates(fragment: Mol, substructure: Mol, atom_matches: Tuple) -> Mol:
    """
    Function to add 3D coordinates to a substructure (e.g. MCS) from the corresponding
    atoms from the original fragment.
    The resulting molecule will be used for constrained embedding.

    :param fragment: the original fragment (with 3D coordinates)
    :type fragment: RDKit molecule
    :param substructure: substructure to add coordinates to
    :type substructure: RDKit molecule
    :param atom_matches: atom numbers of the fragment that match the substructure
    :type atom_matches: tuple

    :return: substructure with coordinates added
    :rtype: RDKit rwmol
    """
    sub_mol = Chem.RWMol(substructure)  # create editable copy of the substructure
    sub_conf = Chem.Conformer(sub_mol.GetNumAtoms())  # conformer of the substructure
    sub_matches = sub_mol.GetSubstructMatch(substructure)  # so atoms in the same order
    ref_conf = fragment.GetConformer()  # get the conformation of the actual fragment

    for i, match in enumerate(
        sub_matches
    ):  # set atom position using matching atom from fragment
        sub_conf.SetAtomPosition(match, ref_conf.GetAtomPosition(atom_matches[i]))

    sub_mol.AddConformer(sub_conf)  # add the conformation to the substructure
    return sub_mol


def remove_xe(synthon: Mol) -> Mol:
    """
    Function to remove the xenon atom from the synthon.

    :param synthon: synthon with xenon denoting attachment point
    :type synthon: RDKit molecule

    :return: synthon with xenon removed
    :rtype: RDKit molecule
    """
    xe = Chem.MolFromSmiles("[Xe]")
    synth = AllChem.DeleteSubstructs(synthon, xe)
    return synth


class EmbeddingFilter(Filter_generic):
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
    def _get_mcs(full_mol: Mol, fragment: Mol) -> Mol:
        """
        Function to return the MCS between a mol and a fragment.

        :param full_mol: full molecule (the merge)
        :type full_mol: RDKit molecule
        :param fragment: fragment the molecule was generated from
        :type fragment: RDKit molecule

        :return: molecule representing the MCS
        :rtype: RDKit molecule (from smarts)
        """
        mcs = rdFMCS.FindMCS([full_mol, fragment], completeRingsOnly=True)
        mcs_mol = Chem.MolFromSmarts(mcs.smartsString)
        Chem.SanitizeMol(mcs_mol)
        return mcs_mol

    @staticmethod
    def get_distance(coord1: np.ndarray, coord2: np.ndarray) -> float:
        """
        Function calculates the distance between two atoms in 3D space.
        Relevant for when two fragments are overlapping.

        :param coord1: 3D atom coordinates
        :type coord1: numpy array
        :param coord2: 3D atom coordinates
        :type coord2: numpy array

        :return: distance between the coordinates
        :rtype: float
        """
        sq = (coord1 - coord2) ** 2
        return np.sqrt(np.sum(sq))

    def _check_overlap(
        self, molA: Mol, molB: Mol, atom_clash: float = config_filter.ATOM_CLASH_DIST
    ) -> Tuple[Mol, Mol]:
        """
        Function checks if parts of two molecules overlap. If atoms overlap, then they are removed
        from one of the molecules.

        :param molA: substructure of first fragment
        :type molA: RDKit molecule
        :param molB: substructure of second fragment
        :type molB: RDKit molecule

        :return: tuple of RDKit molecules with overlapping atoms removed
        :rtype: tuple
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
                dist = self.get_distance(posA, posB)
                # if atoms closer than e.g. 0.5A, recognised as clashes
                if dist < atom_clash:
                    clashes.append(i)
                    break
        if clashes:  # if any atoms identified as clashes
            # sort in descending order (so atoms removed correctly)
            sorted_list = sorted(clashes, reverse=True)
            for clash in sorted_list:
                A.RemoveAtom(clash)  # remove atoms from one fragment
        return A, B

    def embedding(self, fragA: Mol, fragB: Mol, merge_mol: Mol, synth: Mol) -> list:
        """
        Function to embed the full molecule, constraining the atoms that came from each fragment.
        The atoms that came from each fragment are retreived, and the 3D coordinates
        are added from the original structures.

        :param fragA: original fragment A with 3D conformation
        :type fragA: RDKit molecule
        :param fragB: original fragment B with 3D conformation
        :type fragB: RDKit molecule
        :param merge_mol: proposed merge
        :type merge_mol: RDKit molecule
        :param synth: synthon from fragment B used for expansion
        :type synth: RDKit molecule

        :return: list of embedded molecules (if embedding was successful, otherwise empty)
        :rtype: list
        """
        # substructure match for fragment A
        # identify the atoms that came from fragment A
        mcsA = self._get_mcs(merge_mol, fragA)
        # get all possible matches with fragment A
        mcsA_matches = fragA.GetSubstructMatches(mcsA)
        # for each, add coordinates
        mcsA_mols = [add_coordinates(fragA, mcsA, matches) for matches in mcsA_matches]

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
                mcsA_mol, synthB_mol = self._check_overlap(_mcsA_mol, _synthB_mol)
                # combine substructures to get ref molecule
                ref_mol = Chem.CombineMols(mcsA_mol, synthB_mol)
                ref_mols.append(ref_mol)

        # Get substructure matches for merge and embed with all sets of coordinates
        merge_matches = merge_mol.GetSubstructMatches(ref_mols[0])
        embedded_mols = []
        n_embeddings = 0
        for matches in merge_matches:
            for ref_mol in ref_mols:
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

    @staticmethod
    def calc_energy(mol: Mol) -> float:
        """
        Function to calculate the energy of the embedded molecule.

        :param mol: embedded molecule
        :type mol: RDKit molecule

        :return: energy of the molecule
        :rtype: float
        """
        mol_energy = AllChem.UFFGetMoleculeForceField(mol).CalcEnergy()
        return mol_energy

    def calc_unconstrained_energy(
        self, og_mol: Mol, n_conf: int = config_filter.N_CONFORMATIONS
    ) -> float:
        """
        Create ten unconstrained conformations for each molecule and calculate the energy.

        :param og_mol: the original merge without coordinates added
        :type og_mol: RDKit molecule
        :param n_conf: the number of unconstrained conformations to generate
        :type n_conf: int

        :return: the average of the unconstrained energies
        :rtype: float
        """
        unconstrained_energies = []
        for i in range(n_conf):  # generate conformations and calculate energy
            mol = Chem.Mol(og_mol)
            AllChem.EmbedMolecule(mol)
            AllChem.UFFOptimizeMolecule(mol)
            e = self.calc_energy(mol)
            unconstrained_energies.append(e)
        # calculate the average of all the energies
        avg = sum(unconstrained_energies) / len(unconstrained_energies)
        return avg

    def filter_smi(
        self,
        merge: str,
        fragA: Mol,
        fragB: Mol,
        synthon: str,
        energy_threshold: float = config_filter.ENERGY_THRESHOLD,
    ) -> Tuple[bool, None]:
        """
        Runs the filter by embedding (if possible) the molecule, calculating the energy of
        the constrained conformation, and comparing with the energy of the unconstrained
        conformations (by averaging over 10 unconstrained conformations). If the energy ratio is
        greater than a specified threshold for the constrained molecule, the molecule is filtered out.

        :param merge: merge SMILES string
        :type merge: str
        :param fragA: fragment A molecule
        :type fragA: RDKit molecule
        :param fragB: fragment B molecule
        :type fragB: RDKit molecule
        :param synthon: synthon SMILES string
        :type synthon: str
        :param energy_threshold: threshold of constrained energy versus unconstrained energy ratio
        :type energy_threshold: float

        :return: whether molecule passes (True) or fails (False) filter
        :rtype: bool
        """
        merge = Chem.MolFromSmiles(merge)
        synthon = Chem.MolFromSmiles(synthon)
        embedded_mols = self.embedding(fragA, fragB, merge, synthon)

        result, embedded = False, None
        if len(embedded_mols) == 0:
            result = False
            embedded = None

        else:
            for embedded_mol in embedded_mols:
                # energy of constrained conformation
                const_energy = self.calc_energy(embedded_mol)
                # energy of avg unconstrained conformation
                unconst_energy = self.calc_unconstrained_energy(merge)
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

    def filter_all(
        self, cpus: int = config_filter.N_CPUS_FILTER_PAIR
    ) -> Tuple[list, list]:
        """
        Runs the descriptor filter on all the SMILES in parallel.

        :param cpus: number of CPUs for parallelization
        :type cpus: int

        :return: list of results (True or False); list of mols (RDKit molecules)
        :rtype: tuple
        """
        res = Parallel(n_jobs=cpus, backend="multiprocessing")(
            delayed(self.filter_smi)(smi, self._fragmentA, self._fragmentB, synthon)
            for smi, synthon in zip(self.smis, self.synthons)
        )
        self.results = [r[0] for r in res]
        self.mols = [r[1] for r in res]
        return self.results, self.mols
