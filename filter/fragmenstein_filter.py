"""Place and filter the smiles with Fragmenstein"""

import pyrosetta
import json
import os
import shutil
import tempfile

from rdkit import Chem
from rdkit.Chem import rdmolfiles
from joblib import Parallel, delayed
from fragmenstein import Victor
from typing import Tuple

from filter.config_filter import config_filter
from filter.generic_filter import Filter_generic


class FragmensteinFilter(Filter_generic):

    def __init__(self, smis: list, synthons: list, fragmentA, fragmentB, proteinA, proteinB,
                 merge=None, mols=None):
        super().__init__(smis, synthons, fragmentA, fragmentB, proteinA, proteinB, merge, mols)
        self.results = None

    def _get_dict(self, json_file: str) -> dict:
        """
        Function opens the json file to load the dictionary

        :param json_file: json containing the dictionary
        :type json_file: .json
        :return: dictionary containing Fragmenstein info
        :rtype: nested dictionary
        """
        with open(json_file) as f:
            data = json.load(f)
        return data

    def _create_directories(self, working_dir: str = config_filter.WORKING_DIR,
                            output_dir: str = config_filter.OUTPUT_DIR):
        """
        Create subdirectories to store the temporary files and Fragmenstein results.

        :param output_dir: where to save final Fragmenstein files
        :type output_dir: str
        :param working_dir: where to save intermediate tempfiles
        :type working_dir: str
        """
        # create output directory
        if not os.path.exists(output_dir):
            os.mkdir(output_dir)

        # create tempfiles folder in working directory to save intermediate files
        if not os.path.exists((os.path.join(working_dir, 'tempfiles'))):
            os.mkdir(os.path.join(working_dir, 'tempfiles'))

        # create Fragmenstein folder to save the important files we want to keep
        if not os.path.exists(os.path.join(output_dir, 'fragmenstein')):
            os.mkdir(os.path.join(output_dir, 'fragmenstein'))

    def _get_name(self, num: int) -> str:
        """
        Name the files using the names of the fragments (get from their files), e.g. x0034_0B_x0176_0B_24

        :param num: identified to use for naming the merge
        :type num: int

        :return: name
        :rtype: str
        """
        name = self.merge.replace('_', '-')
        name = f'{name}-{num}'
        return name

    def _move_file(self, name, last_dir, new_dir, ext):
        # move the json file
        new_file = name + ext
        old_path = os.path.join(last_dir, name, new_file)
        new_path = os.path.join(new_dir, 'fragmenstein', new_file)
        shutil.move(old_path, new_path)
        return new_path

    def place_smiles(self, name, smiles, working_dir: str = config_filter.WORKING_DIR,
                     output_dir: str = config_filter.OUTPUT_DIR, residue=config_filter.COVALENT_RESI):
        """

        Function places the merge in the protein using Fragmenstein. Minimization is performed
        using PyRosetta. Creates an output folder containing .json file with metrics that
        can be used for filtering.

        :param name: name of the merge (to name the output subfolder/files)
        :type name: str
        :param smiles: smiles of the merge
        :type smiles: str
        :param working_dir: filepath of working directory
        :type working_dir: str
        :param output_dir: filepath of output directory
        :type output_dir: str
        :param residue: covalent residue for PyRosetta, e.g. 2B
        :type residue: str
        """
        # create temporary directory to write files to
        temp_dir = tempfile.TemporaryDirectory(dir=f'{working_dir}/tempfiles/')

        # run Fragmenstein
        # initialise PyRosetta
        pyrosetta.init(
            extra_options='-no_optH false -mute all -ex1 -ex2 -ignore_unrecognized_res false -load_PDB_components false -ignore_waters false')

        # get the fragments from the files
        fragments_fnames = [self.fragmentA, self.fragmentB]
        hits = [Chem.MolFromMolFile(frag) for frag in fragments_fnames]

        # set the output directory
        Victor.work_path = temp_dir.name

        # set up Victor and place the smiles
        v = Victor(hits=hits, pdb_filename=self.proteinA, covalent_resi=residue)
        v.place(smiles=smiles, long_name=name)  # 'long_name' used to name the files

        save_files = ['.minimised.json', '.minimised.mol', '.holo_minimised.pdb']
        new_files = []
        for f in save_files:
            new_f = self._move_file(name, temp_dir.name, output_dir, f)
            new_files.append(new_f)

        # remove files from temporary directory
        temp_dir.cleanup()

        return new_files

    def filter_smi(self, idx: int, smi: str, comRMSD_threshold: float = config_filter.COM_RMSD) -> Tuple:
        """
        Function filters molecules for those where both fragments were considered in its placement,
        a negative ΔΔG and combined RMSD with the fragments of < 1.5A.

        :param idx: number of SMILES for naming the file
        :type idx: int
        :param smi: merge SMILES string
        :type smi: str
        :param comRMSD_threshold: threshold for combined RMSD filter
        :type comRMSD_threshold: float

        :return: returns 'pass' or 'fail'; molecule minimized by Fragmenstein
        :rtype: tuple
        """
        merge_name = self._get_name(idx)
        new_files = self.place_smiles(merge_name, smi)

        # get json_file for filtering and mol_file to get mol
        json_file = new_files[0]
        mol_file = new_files[1]
        minimised_mol = rdmolfiles.MolFromMolFile(mol_file)

        data = self._get_dict(json_file)  # load the dictionary from the json file

        # retrieve the energy of the bound and unbound molecules and the comRMSD
        G_bound = data['Energy']['ligand_ref2015']['total_score']  # energy of bound molecule
        G_unbound = data['Energy']['unbound_ref2015']['total_score']  # energy of unbound molecule
        deltaG = G_bound - G_unbound  # calculate energy difference
        comRMSD = data['mRMSD']  # RMSD between two fragments and merge

        # get number of fragments used for placement of SMILES
        regarded = 0
        for rmsd in data['RMSDs']:
            if rmsd:
                regarded += 1

        # only keep molecules where both fragments used in placement, -ΔΔG and comRMSD < threshold
        if (regarded == 2) & (deltaG < 0) & (comRMSD <= comRMSD_threshold):
            result = True
        else:
            result = False
            for f in new_files:
                os.remove(f)

        return result, minimised_mol

    def filter_all(self, cpus: int = config_filter.N_CPUS_FILTER_PAIR) -> Tuple[list, list]:
        """
        Runs the Fragmenstein filter on all the SMILES in parallel.

        :param cpus: number of CPUs for parallelization
        :type cpus: int

        :return: list of results (True or False); list of mols (RDKit molecules)
        :rtype: tuple
        """
        self._create_directories()
        idxs = range(len(self.smis))  # generate numbers to number results in their filenames
        res = Parallel(n_jobs=cpus, backend='multiprocessing')(delayed(self.filter_smi)(idx, smi) for idx, smi in
                                                               zip(idxs, self.smis))
        self.results = [r[0] for r in res]
        self.mols = [r[1] for r in res]
        return self.results, self.mols
