"""
Abstract class for scoring step in the pipeline
"""

import os
import shutil
from abc import ABC, abstractmethod

from fragment_network_merges.filter.config_filter import config_filter
from joblib import Parallel, delayed
from rdkit.Chem import rdmolfiles
from fragment_network_merges.utils.filter_utils import remove_ligand


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

    def get_apo_files(self, cpus: int = config_filter.N_CPUS_FILTER_PAIR):
        """
        Get the apo files by removing ligand from holo files.
        """
        self.apo_files = Parallel(n_jobs=cpus, backend="multiprocessing")(
            delayed(remove_ligand)(holo_file) for holo_file in self.holo_files
        )

    def move_apo_files(self):
        """
        Move the apo files to the output directory
        """
        for _name, pdb_file in zip(self.names, self.apo_files):
            name = _name.replace("_", "-")
            fname = os.path.basename(pdb_file)
            new_path = os.path.join(self.out_pair_dir, name, fname)
            shutil.copy(pdb_file, new_path)

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
