"""
Abstract class for filtering step in the pipeline
"""

from abc import ABC, abstractmethod

from rdkit.Chem import rdmolfiles


class Filter_generic(ABC):
    """
    Abstract class for filtering step
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
    ):
        self.names = names
        self.smis = smis  # list of SMILES of merges
        self.synthons = synthons  # list of synthons corresponding to SMILES of merges
        self.fragmentA = fragmentA  # filepath
        self.fragmentB = fragmentB  # filepath
        self.proteinA = proteinA  # filepath
        self.proteinB = proteinB  # filepath
        self.merge = merge  # SMILES representing the merge (e.g. x0001_0B_x0002_0B)
        self.mols = mols  # list of molecules with conformers (if generated)
        self.work_pair_dir = work_pair_dir
        self.out_pair_dir = out_pair_dir

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
        return self.mol_files, self.holo_files, self.apo_files

    @abstractmethod
    def filter_smi(self, **kwargs):
        """
        Passes a single SMILES through the filter.
        """
        raise NotImplementedError

    @abstractmethod
    def filter_all(self, **kwargs):
        """
        Runs the filter on all SMILES for a merge pair in parallel.
        """
        raise NotImplementedError
