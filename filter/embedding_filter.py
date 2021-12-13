"""
Used for 3D filtering of fragment merges by constrained embedding.
"""

import numpy as np

from joblib import Parallel, delayed
from rdkit import Chem
from rdkit import RDLogger
from rdkit.Chem import AllChem, rdFMCS, rdForceFieldHelpers
from typing import Tuple

from filter.config_filter import config_filter
from filter.generic_filter import Filter_generic

RDLogger.DisableLog('rdApp.*')


def add_coordinates(fragment, substructure):
    """
    Function to add 3D coordinates to a substructure (e.g. MCS) from the corresponding
    atoms from the original fragment.
    The resulting molecule will be used for constrained embedding.

    :param fragment: the original fragment (with 3D coordinates)
    :type fragment: RDKit molecule
    :param substructure: substructure to add coordinates to
    :type substructure: RDKit molecule

    :return: substructure with coordinates added
    :rtype: RDKit rwmol
    """
    ref_match = fragment.GetSubstructMatch(substructure)  # get frag atoms that match substruct
    rwmol = Chem.RWMol(substructure)  # create editable copy of the substructure
    rwconf = Chem.Conformer(rwmol.GetNumAtoms())  # create a conformer of the substructure
    matches = rwmol.GetSubstructMatch(substructure)  # get matches so atoms in the same order
    ref_conf = fragment.GetConformer()  # get the conformation of the actual fragment
    for i, match in enumerate(matches):  # set atom position using matching atom from fragment
        rwconf.SetAtomPosition(match, ref_conf.GetAtomPosition(ref_match[i]))
    rwmol.AddConformer(rwconf)  # add the conformation to the substructure
    return rwmol


def remove_xe(synthon):
    """
    Function to remove the xenon atom from the synthon.

    :param synthon: synthon with xenon denoting attachment point
    :type synthon: RDKit molecule

    :return: synthon with xenon removed
    :rtype: RDKit molecule
    """
    xe = Chem.MolFromSmiles('[Xe]')
    synth = AllChem.DeleteSubstructs(synthon, xe)
    return synth


class EmbeddingFilter(Filter_generic):

    def __init__(self, smis: list, synthons: list, fragmentA, fragmentB, proteinA=None, proteinB=None, merge=None,
                 mols=None):
        super().__init__(smis, synthons, fragmentA, fragmentB, proteinA, proteinB, merge, mols)
        self.results = None

    def _get_mcs(self, full_mol, fragment):
        """
        Function to return the MCS between a mol and a fragment.

        :param full_mol: full molecule (the merge)
        :type full_mol: RDKit molecule
        :param fragment: fragment the molecule was generated from
        :type fragment: RDKit molecule

        :return: molecule representing the MCS
        :rtype: RDKit molecule (from smarts)
        """
        mcs = rdFMCS.FindMCS([full_mol, fragment])
        mcs_mol = Chem.MolFromSmarts(mcs.smartsString)
        return mcs_mol

    def get_distance(self, coord1: np.ndarray, coord2: np.ndarray) -> float:
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

    def _check_overlap(self, molA, molB):
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
                # if atoms closer than 0.5A, recognised as clashes
                if dist < 0.5:
                    clashes.append(i)
                    break
        if clashes:  # if any atoms identified as clashes
            # sort in descending order (so atoms removed correctly)
            sorted_list = sorted(clashes, reverse=True)
            for clash in sorted_list:
                A.RemoveAtom(clash)  # remove atoms from one fragment
        return A, B

    def embedding(self, fragA, fragB, full_mol, synth):
        """
        Function to embed the full molecule, constraining the atoms that came from each fragment.
        The atoms that came from each fragment are retreived, and the 3D coordinates
        are added from the original structures.

        :param fragA: original fragment A with 3D conformation
        :type fragA: RDKit molecule
        :param fragB: original fragment B with 3D conformation
        :type fragB: RDKit molecule
        :param full_mol: proposed merge
        :type full_mol: RDKit molecule
        :param synth: synthon from fragment B used for expansion
        :type synth: RDKit molecule

        :return: embedded molecule (if embedding was successful)
        :rtype: RDKit molecule
        """
        mcsA = self._get_mcs(full_mol, fragA)  # identify the atoms that came from fragment A
        synthB = remove_xe(synth)  # remove the xenon from the synthon
        rwmolA = add_coordinates(fragA, mcsA)
        rwmolB = add_coordinates(fragB, synthB)
        newmolA, newmolB = self._check_overlap(rwmolA, rwmolB)  # check if atoms overlap
        combined_mol = Chem.CombineMols(newmolA, newmolB)  # combine mols to get reference molecule
        # full_mol = Chem.AddHs(full_mol)  #TODO: test with H removal
        embedded = AllChem.ConstrainedEmbed(Chem.Mol(full_mol), combined_mol, 42)  # do embedding
        rdForceFieldHelpers.MMFFOptimizeMolecule(embedded)  # optimize the embedding
        # embedded = Chem.RemoveHs(embedded)
        return embedded

    def calc_energy(self, mol) -> float:
        """
        Function to calculate the energy of the embedded molecule.

        :param mol: embedded molecule
        :type mol: RDKit molecule

        :return: energy of the molecule
        :rtype: float
        """
        mol_energy = AllChem.UFFGetMoleculeForceField(mol).CalcEnergy()
        return mol_energy

    def calc_unconstrained_energy(self, og_mol) -> float:
        """
        Create ten unconstrained conformations for each molecule and calculate the energy.

        :param og_mol: the original merge without coordinates added
        :type og_mol: RDKit molecule

        :return: the average of the unconstrained energies
        :rtype: float
        """
        unconstrained_energies = []
        for i in range(10):  # generate 10 conformations and calculate energy
            mol = Chem.Mol(og_mol)
            AllChem.EmbedMolecule(mol)
            AllChem.UFFOptimizeMolecule(mol)
            e = self.calc_energy(mol)
            unconstrained_energies.append(e)
        # calculate the average of all the energies
        avg = sum(unconstrained_energies) / len(unconstrained_energies)
        return avg

    def filter_smi(self, merge: str, fragA, fragB, synthon: str, energy_threshold: float = config_filter.ENERGY_THRESHOLD)\
            -> Tuple:
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
        try:
            embedded = self.embedding(fragA, fragB, merge, synthon)
            const_energy = self.calc_energy(embedded)
            unconst_energy = self.calc_unconstrained_energy(merge)
            # if the energy of the constrained conformation is less, then pass filter
            if const_energy <= unconst_energy:
                result = True
            else:
                # if constrained energy >energy-threshold-fold greater, then fail filter
                ratio = const_energy / unconst_energy
                if ratio >= energy_threshold:
                    result = False
                    embedded = None
                else:
                    result = True
        except:
            result = False  # if embedding fails, then fail filter
            embedded = None

        return result, embedded

    def filter_all(self, cpus: int = config_filter.N_CPUS_FILTER_PAIR) -> Tuple[list, list]:
        """
        Runs the descriptor filter on all the SMILES in parallel.

        :param cpus: number of CPUs for parallelization
        :type cpus: int

        :return: list of results (True or False); list of mols (RDKit molecules)
        :rtype: tuple
        """
        res = Parallel(n_jobs=cpus, backend='multiprocessing') \
            (delayed(self.filter_smi)(smi, self._fragmentA, self._fragmentB, synthon) for smi, synthon in
             zip(self.smis, self.synthons))
        self.results = [r[0] for r in res]
        self.mols = [r[1] for r in res]
        return self.results, self.mols
