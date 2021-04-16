"""
Used for 3D filtering of fragment merges by constrained embedding.
"""

import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem, rdFMCS, rdForceFieldHelpers

class EmbeddingFilter():
    """Filters molecules using constrained embedding using coordinates from original fragments."""

    def __init__(self, merge, fragA, fragB, synthon):
        self.merge = Chem.MolFromSmiles(merge)
        self.fragA = fragA
        self.fragB = fragB
        self.synthon = Chem.MolFromSmiles(synthon)
        self.embedded = None
        self.result = None

    def get_mcs(self, full_mol, fragment):
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

    def add_coordinates(self, fragment, substructure):
        """
        Function to add 3D coordinates to a substructure (e.g. MCS) from the corresponding
        atoms from the original fragment.
        The resulting molecule will be used for constrained embedding.

        :param fragment: the original fragment with 3D coordinates
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
            # Added atom position information from reference molecule
            rwconf.SetAtomPosition(match, ref_conf.GetAtomPosition(ref_match[i]))
        rwmol.AddConformer(rwconf)  # add the conformation to the substructure
        return rwmol

    def get_distance(self, coord1, coord2):
        """
        Function calculates the distance between two atoms in 3D space.
        Relevant for when two fragments are overlapping.
        Distance calculated using Pythagoras.

        :param coord1: atom coordinates
        :type coord1: 3D coordinates
        :param coord2: atom coordinates
        :type coord2: 3D coordinates

        :return: distance between the coordinates
        :rtype: float
        """
        sq = (coord1 - coord2) ** 2
        return np.sqrt(np.sum(sq))

    def check_overlap(self, molA, molB):
        """
        Function checks if parts of two molecules overlap. If atoms overlap, then they are removed
        from one of the molecules.

        :param molA: substructure of first fragment
        :type molA: RDKit molecule
        :param molB: substructure of second fragment
        :type molB: RDKit molecule

        :return: two molecules with overlapping atoms removed
        :rtype: RDKit molecules
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
                if dist < 0.5:
                    clashes.append(i)
                    break
        if clashes:
            s = sorted(clashes, reverse=True)
            for c in s:
                A.RemoveAtom(c)
        return A, B

    def remove_xe(self, synthon):
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
        mcsA = self.get_mcs(full_mol, fragA)
        synthB = self.remove_xe(synth)
        rwmolA = self.add_coordinates(fragA, mcsA)
        rwmolB = self.add_coordinates(fragB, synthB)
        newmolA, newmolB = self.check_overlap(rwmolA, rwmolB) # check if atoms overlap
        combined_mol = Chem.CombineMols(newmolA, newmolB) # combine mols to get reference molecule
        embedded = AllChem.ConstrainedEmbed(Chem.Mol(full_mol), combined_mol, 42) # do embedding
        rdForceFieldHelpers.MMFFOptimizeMolecule(embedded) # optimize the embedding
        return embedded

    def calc_energy(self, mol):
        """
        Funcion to calculate the energy of the embedded molecule.

        :param mol: embedded molecule
        :type mol: RDKit molecule

        :return: energy of the molecule
        :rtype: float
        """
        mol_energy = AllChem.UFFGetMoleculeForceField(mol).CalcEnergy()
        return mol_energy

    def calc_unconstrained_energy(self, og_mol):
        """
        Create ten unconstrained conformations for each molecule and calculate the energy.

        :param og_mol: the original merge without coordinates added
        :type og_mol: rdkit molecule

        :return: the average of the unconstrained energies
        :rtype: float
        """
        unconstrained_energies = []
        for i in range(10):
            mol = Chem.Mol(og_mol)
            AllChem.EmbedMolecule(mol)
            AllChem.UFFOptimizeMolecule(mol)
            e = self.calc_energy(mol)
            unconstrained_energies.append(e)
        avg = sum(unconstrained_energies) / len(unconstrained_energies)
        return avg

    def filter(self):
        """
        Runs the filter by embedding (if possible) the molecule, calculating the energy of
        the constrained conformation, and comparing with the energy of the unconstrained
        conformations (by averaging over 10 unconstrained conformations). If the energy is
        >10-fold greater for the constrained molecule, the molecule is filtered out.

        :return: embedded molecule or None
        :rtype: rdkit molecule or None
        """
        try:
            self.embedded = self.embedding(self.fragA, self.fragB, self.merge, self.synthon)
            const_energy = self.calc_energy(self.embedded)
            unconst_energy = self.calc_unconstrained_energy(self.merge)
            if const_energy <= unconst_energy:
                self.result = 'pass'
            else:
                ratio = const_energy / unconst_energy
                if ratio >= 10:
                    self.result = 'fail'
                else:
                    self.result = 'pass'
        except:
            self.result = 'fail'

        return self.result
