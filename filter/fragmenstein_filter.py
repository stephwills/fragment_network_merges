"""Place and filter the smiles with Fragmenstein"""

import json
import multiprocessing as mp
import multiprocessing.queues as mpq
from multiprocessing.queues import Queue
import os
import shutil
import time
from typing import Tuple

import pyrosetta
from filter.config_filter import config_filter
from filter.generic_filter import Filter_generic
from fragmenstein import Victor
from rdkit import Chem
from rdkit.Chem import rdmolfiles


class FragmensteinFilter(Filter_generic):
    def __init__(
        self,
        smis: list,
        synthons: list,
        fragmentA,
        fragmentB,
        proteinA,
        proteinB,
        merge,
        mols=None,
        names=None,
        work_pair_dir=None,
        out_pair_dir=None,
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
        )
        self.results = None
        self.timings = None

    def _copy_files(self, name: str, merge_dir: str):
        """
        If filter passes, move the working dir files to the output dir.

        :param name: name of the merge used in file naming
        :type name: str
        :param merge_dir: name of the merge directory
        :type merge_dir: str
        """
        output_subdir = os.path.join(self.out_pair_dir, name)
        if not os.path.exists(output_subdir):
            os.mkdir(output_subdir)
        dir_files = os.listdir(merge_dir)
        workdir_files = [os.path.join(merge_dir, f) for f in dir_files]
        outputdir_files = [os.path.join(output_subdir, f) for f in dir_files]
        for workdir_file, outputdir_file in zip(workdir_files, outputdir_files):
            shutil.copy(workdir_file, outputdir_file)

    @staticmethod
    def _get_dict(json_file: str) -> dict:
        """
        Function opens the json file to load the dictionary

        :param json_file: json containing the dictionary
        :type json_file: .json

        :return: dictionary containing Fragmenstein info
        :rtype: nested dictionary
        """
        if os.path.exists(json_file):
            with open(json_file) as f:
                data = json.load(f)
            return data

    def _check_merge_directory(self, name: str) -> Tuple[str, bool]:
        """
        Create sub-directories for each merge to save the fragmenstein files to. These will all be saved in the
        working directory and then moved later if filtering is successful (to avoid reading and writing to shared
        disk).

        :param name: unique name of the merge
        :type name: str

        :return: merge_dir, bool (True if fragmenstein not run, False if already run)
        :rtype: tuple
        """
        merge_dir = os.path.join(self.work_pair_dir, name)
        fname = f"{name}.minimised.json"
        fname = os.path.join(merge_dir, fname)
        if os.path.exists(fname):
            return merge_dir, True
        else:
            return merge_dir, False

    @staticmethod
    def _get_name(name: str) -> str:
        """
        Get the merge names used in the file paths - Fragmenstein replaces the underscores with hyphens.

        :return: name
        :rtype: str
        """
        _name = name.replace("_", "-")
        return _name

    def filter_smi(
        self,
        queue: Queue,
        merge_name: str,
        smi: str,
        comRMSD_threshold: float = config_filter.COM_RMSD,
        residue: int = config_filter.COVALENT_RESI,
    ):
        """
        Function filters molecules for those where both fragments were considered in its placement,
        a negative ΔΔG and combined RMSD with the fragments of < 1.5A. queue.put returns list with whether the molecule
        passes and the minimised mol, holo and apo files if successful.

        :param queue: multiprocessing queue to allow timeout
        :type queue: multiprocessing queue
        :param merge_name: name of the unique merge (e.g. x0107_0A_x0034_0A_24)
        :type merge_name: str
        :param smi: merge SMILES string
        :type smi: str
        :param comRMSD_threshold: threshold for combined RMSD filter
        :type comRMSD_threshold: float
        :param residue: covalent residue for PyRosetta, e.g. 2B
        :type residue: str
        """
        start = time.time()

        name = self._get_name(merge_name)
        idx = int(name.split('-')[-1])
        # check if directory already exists within the working directory to write intermediate files to
        merge_dir, already_run = self._check_merge_directory(name)

        json_file = os.path.join(merge_dir, f"{name}.minimised.json")

        if not already_run:
            # run Fragmenstein
            # initialise PyRosetta
            pyrosetta.init(
                extra_options="""-no_optH false -mute all -ex1 -ex2 -ignore_unrecognized_res false
                                            -load_PDB_components false -ignore_waters false"""
            )

            # get the fragments from the files
            fragments_fnames = [self.fragmentA, self.fragmentB]
            hits = [Chem.MolFromMolFile(frag) for frag in fragments_fnames]

            # set the output directory
            Victor.work_path = self.work_pair_dir

            try:
                # set up Victor and place the smiles
                v = Victor(hits=hits, pdb_filename=self.proteinA, covalent_resi=residue)
                v.place(smiles=smi, long_name=name)  # 'long_name' used to name the files
            except Exception as e:
                # if Fragmenstein fails, record the error and return filter fail
                print(f"Fragmenstein failed for {name}: {str(e)}")
                self.errors[name] = str(e)
                with open(
                        os.path.join(
                            self.work_pair_dir, f"{self.merge}_fragmenstein_errors.json"
                        ),
                        "w",
                ) as f:
                    json.dump(self.errors.copy(), f)
                queue.put([idx, False, None, None, None])

        data = self._get_dict(json_file)

        if data:
            # retrieve the energy of the bound and unbound molecules and the comRMSD
            G_bound = data["Energy"]["ligand_ref2015"][
                "total_score"
            ]  # energy of bound molecule
            G_unbound = data["Energy"]["unbound_ref2015"][
                "total_score"
            ]  # energy of unbound molecule
            deltaG = G_bound - G_unbound  # calculate energy difference
            comRMSD = data["mRMSD"]  # RMSD between two fragments and merge

            # get number of fragments used for placement of SMILES
            regarded = 0
            for rmsd in data["RMSDs"]:
                if rmsd:
                    regarded += 1

            # get other files
            mol_file = os.path.join(merge_dir, f"{name}.minimised.mol")
            holo_file = os.path.join(merge_dir, f"{name}.holo_minimised.pdb")
            minimised_mol = rdmolfiles.MolFromMolFile(mol_file)

            # only keep molecules where both fragments used in placement, -ΔΔG and comRMSD < threshold
            if (regarded == 2) & (deltaG < 0) & (comRMSD <= comRMSD_threshold):
                result = True
                # move the fragmenstein files to the output directory
                self._copy_files(name, merge_dir)

            else:
                result = False
                # if filter failed, remove all those files from the working directory folder
                shutil.rmtree(merge_dir)

            # record timings for fragmenstein
            end = time.time()
            total = round(end - start, 2)
            timings_dict = {'timings': total, 'result': result}
            self.timings[name] = timings_dict
            with open(
                os.path.join(
                    self.work_pair_dir, f"{self.merge}_fragmenstein_timings.json"
                ),
                "w",
            ) as f:
                json.dump(self.timings.copy(), f)
            queue.put([idx, result, minimised_mol, mol_file, holo_file])

        else:
            queue.put([idx, False, None, None, None])

    def filter_all(
        self, cpus: int = config_filter.N_CPUS_FILTER_PAIR
    ) -> Tuple[list, list]:
        """
        Runs the Fragmenstein filter on all the SMILES in parallel.

        :param cpus: number of CPUs for parallelization
        :type cpus: int

        :return: list of results (True or False); list of mols (RDKit molecules)
        :rtype: tuple
        """
        # create multiprocessing Process to implement timeout

        queue = mp.Queue()
        manager = mp.Manager()
        self.timings = manager.dict()
        self.errors = manager.dict()
        processes = [
            mp.Process(target=self.filter_smi, args=(queue, name, smi))
            for name, smi in zip(self.names, self.smis)
        ]

        res = []
        for p in processes:
            p.start()

        for p in processes:
            p.join()

        for p in processes:
            try:
                r = queue.get(timeout=5)  # timeout set to 10 minutes
                res.append(r)
            except mpq.Empty:
                p.terminate()
        print('res', res)
        res = sorted(res, key=lambda tup: tup[0])
        self.results = [r[1] for r in res]
        self.mols = [r[2] for r in res]
        self.mol_files = [r[3] for r in res]
        self.holo_files = [r[4] for r in res]

        return self.results, self.mols
