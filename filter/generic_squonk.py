"""
Abstract class for using informatics matters data manager job utilities
"""

import importlib
import time
from abc import ABC

try:
    from dm_job_utilities.dm_log import DmLog
except ImportError:
    class DmLog():
        @staticmethod
        def emit_event(*args, **kwargs):
            pass
        @staticmethod
        def emit_cost(*args, **kwargs):
            pass
from filter.config_filter import config_filter
from rdkit import Chem


class Squonk_generic(ABC):
    """
    Abstract class for running filter as squonk job
    """

    def __init__(
        self,
        filter_name,
        args,
        input_sdf,
        output_sdf,
        fragmentA=None,
        fragmentB=None,
        proteinA=None,
        proteinB=None,
        work_pair_dir=None,
        out_pair_dir=None,
        pair=None,
        score=None,
    ):
        self.filter_name = filter_name
        self.args = args
        self.input_sdf = input_sdf
        self.output_sdf = output_sdf
        self.fragmentA = fragmentA
        self.fragmentB = fragmentB
        self.proteinA = proteinA
        self.proteinB = proteinB
        self.work_pair_dir = work_pair_dir
        self.out_pair_dir = out_pair_dir
        self.pair = pair
        self.score = score

    def execute_job(self, *args):
        DmLog.emit_event(self.filter_name, self.args)
        start = time.time()
        count = 0
        hits = 0
        errors = 0

        w = Chem.SDWriter(self.output_sdf)
        suppl = Chem.SDMolSupplier(self.input_sdf)

        mols, smiles, synthons, names, mol_files, holo_files, apo_files = (
            [],
            [],
            [],
            [],
            [],
            [],
            [],
        )
        for mol in suppl:
            if mol is None:
                continue
            else:
                mols.append(mol)
                count += 1
                synthon = None
                mol_file = None
                holo_file = None
                apo_file = None

                if mol.HasProp("name") == 1:
                    name = mol.GetProp("name")
                    names.append(name)
                else:
                    if self.pair:
                        names.append(self.pair + "_" + str(count - 1))

                smiles.append(Chem.MolToSmiles(mol))

                if mol.HasProp("synthon") == 1:
                    synthon = mol.GetProp("synthon")

                if self.score:
                    if mol.HasProp("mol_file") == 1:
                        mol_file = mol.GetProp("mol_file")
                    if mol.HasProp("holo_file") == 1:
                        holo_file = mol.GetProp("holo_file")
                    if mol.HasProp("apo_file") == 1:
                        apo_file = mol.GetProp("apo_file")

                synthons.append(synthon)
                mol_files.append(mol_file)
                holo_files.append(holo_file)
                apo_files.append(apo_file)

        module = config_filter.PIPELINE_DICT[self.filter_name]
        cls = getattr(importlib.import_module(module), self.filter_name)
        filter = cls(
            smis=smiles,
            synthons=synthons,
            fragmentA=self.fragmentA,
            fragmentB=self.fragmentB,
            proteinA=self.proteinA,
            proteinB=self.proteinB,
            mols=mols,
            work_pair_dir=self.work_pair_dir,
            out_pair_dir=self.out_pair_dir,
            names=names,
            merge=self.pair,
        )
        if self.score:
            filter.mol_files = mol_files
            filter.holo_files = holo_files
            filter.apo_files = apo_files

        results, filtered_mols = filter.filter_all(*args)

        for res, filtered_mol, synthon, mol_file, holo_file, apo_file in zip(
            results, filtered_mols, synthons, mol_files, holo_files, apo_files
        ):
            if res:
                hits += 1
                if synthon:
                    filtered_mol.SetProp("synthon", synthon)
                if self.score:
                    if holo_file:
                        filtered_mol.SetProp("holo_file", holo_file)
                    if apo_file:
                        filtered_mol.SetProp("apo_file", apo_file)
                    if mol_file:
                        filtered_mol.SetProp("mol_file", mol_file)
                w.write(filtered_mol)

        w.close()

        end = time.time()
        duration_s = int(end - start)
        if duration_s < 1:
            duration_s = 1
        DmLog.emit_event(
            count, "inputs,", hits, "hits,", errors, "errors.", "Time (s):", duration_s
        )
        DmLog.emit_cost(count)
