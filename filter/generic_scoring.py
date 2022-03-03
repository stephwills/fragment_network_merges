"""
Abstract class for scoring step in the pipeline
"""

import os
import shutil
from abc import ABC, abstractmethod

from Bio.PDB import Select
from Bio.PDB.PDBIO import PDBIO
from Bio.PDB.PDBParser import PDBParser
from rdkit.Chem import rdmolfiles


class ResSelect(Select):
    """
    Class to remove ligand 'residue' from pdb structure.
    """
    def accept_residue(self, residue):
        if residue.id[0] == "H_LIG":
            return False
        else:
            return True


class Score_generic(ABC):
    """
    Abstract class for scoring filtered molecules
    """

    def __init__(
        self,
        smis: list,
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
        self.smis = smis  # list of SMILES of merges
        self.synthons = synthons  # list of synthons corresponding to SMILES of merges
        self.fragmentA = fragmentA  # filepath
        self.fragmentB = fragmentB  # filepath
        self.proteinA = proteinA  # filepath
        self.proteinB = proteinB  # filepath
        self.merge = merge  # SMILES representing the merge (e.g. x0001_0B_x0002_0B)
        self.mols = mols  # list of molecules with conformers (if generated)
        self.names = names  # list of unique merge names (e.g. x0034_0B_x0176_0B_123)
        self.work_pair_dir = work_pair_dir
        self.out_pair_dir = out_pair_dir

        # placed that may be used specifically for scoring
        self.mol_files = mol_files  # list of placed mol files
        self.holo_files = holo_files  # list of placed holo_files
        self.apo_files = apo_files  # list of placed apo_files

        # get mols from filepaths
        self._fragmentA = rdmolfiles.MolFromMolFile(self.fragmentA)  # RDKit molecule
        self._fragmentB = rdmolfiles.MolFromMolFile(self.fragmentB)  # RDKit molecule
        self._proteinA = rdmolfiles.MolFromPDBFile(self.proteinA)  # RDKit molecule
        self._proteinB = rdmolfiles.MolFromPDBFile(self.proteinB)  # RDKit molecule

    def setattrs(self, **kwargs):
        for k, v in kwargs.items():
            setattr(self, k, v)

    def _remove_ligand(self):
        """
        Function uses BioPython to read in pdb structure, remove the ligand
        and save to a new file in the same directory.
        """
        self.apo_files = []
        for pdb_file in self.holo_files:
            parser = PDBParser(PERMISSIVE=1)
            structure = parser.get_structure("pdb", pdb_file)
            new_filename = pdb_file.replace(".pdb", "_nolig.pdb")

            io = PDBIO()
            io.set_structure(structure)

            for model in structure:
                for chain in model:
                    for residue in chain:
                        io.save(new_filename, ResSelect())
            self.apo_files.append(new_filename)

    def move_apo_files(self):
        """
        Move the apo files to the output directory
        """
        for _name, pdb_file in zip(self.names, self.apo_files):
            name = _name.replace('_', '-')
            fname = os.path.basename(pdb_file)
            new_path = os.path.join(self.out_pair_dir, name, fname)
            shutil.move(pdb_file, new_path)

    @abstractmethod
    def score_mol(self, **kwargs):
        """
        Scoring a single filtered molecule.
        """
        raise NotImplementedError

    @abstractmethod
    def score_all(self, **kwargs):
        """
        Scoring all molecules in parallel.
        """
        raise NotImplementedError