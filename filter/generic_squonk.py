"""
Abstract class for using informatics matters data manager job utilities
"""

import importlib
import time
from abc import ABC

from dm_job_utilities.dm_log import DmLog
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
        pair=None
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

    def execute_job(self, *args):
        DmLog.emit_event(self.filter_name, self.args)
        print(self.output_sdf)
        print(self.input_sdf)
        start = time.time()
        count = 0
        hits = 0
        errors = 0

        w = Chem.SDWriter(self.output_sdf)
        suppl = Chem.SDMolSupplier(self.input_sdf)

        mols, smiles, synthons, names = [], [], [], []
        for mol in suppl:
            if mol is None:
                continue
            else:
                count += 1
                names.append(self.pair + '_' + str(count-1))
                smiles.append(Chem.MolToSmiles(mol))
                if mol.HasProp("synthon") == 1:
                    print('synthon', mol.GetProp("synthon"))
                    synthon = mol.GetProp("synthon")
                else:
                    synthon = None
                synthons.append(synthon)
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
            merge=self.pair
        )

        results, filtered_mols = filter.filter_all(*args)

        for res, filtered_mol, synthon in zip(results, filtered_mols, synthons):
            if res:
                hits += 1
                filtered_mol.SetProp("synthon", synthon)
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
