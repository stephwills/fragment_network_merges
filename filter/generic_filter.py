"""
Abstract class for filtering step in the pipeline
"""

from rdkit.Chem import rdmolfiles
from abc import ABC, abstractmethod


class Filter_generic(ABC):
    """
    Abstract class for filtering step
    """
    def __init__(self, smis: list, synthons=None, fragmentA=None, fragmentB=None, proteinA=None,
                 proteinB=None, merge=None, mols=None, names=None):
        self.names = names
        self.smis = smis  # list of SMILES of merges
        self.synthons = synthons  # list of synthons corresponding to SMILES of merges
        self.fragmentA = fragmentA  # filepath
        self.fragmentB = fragmentB  # filepath
        self.proteinA = proteinA  # filepath
        self.proteinB = proteinB  # filepath
        self.merge = merge # SMILES representing the merge (e.g. x0001_0B_x0002_0B)
        self.mols = mols  # list of molecules with conformers (if generated)

        # get mols from filepaths
        self._fragmentA = rdmolfiles.MolFromMolFile(self.fragmentA)  # RDKit molecule
        self._fragmentB = rdmolfiles.MolFromMolFile(self.fragmentB)  # RDKit molecule
        self._proteinA = rdmolfiles.MolFromPDBFile(self.proteinA)  # RDKit molecule
        self._proteinB = rdmolfiles.MolFromPDBFile(self.proteinB)  # RDKit molecule

        # to store filepaths of placed files
        self.mol_files = None
        self.apo_files = None
        self.holo_files = None

    def setattrs(self, **kwargs):
        for k, v in kwargs.items():
            setattr(self, k, v)

    def get_placed_files(self):
        return self.mol_files, self.apo_files, self.holo_files

    @abstractmethod
    def filter_smi(self):
        """
        Passes a single SMILES through the filter.
        """
        raise NotImplementedError

    @abstractmethod
    def filter_all(self):
        """
        Runs the filter on all SMILES for a merge pair in parallel.
        """
        raise NotImplementedError
