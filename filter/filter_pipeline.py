"""Run the filtering pipeline"""
import argparse
import importlib
import os
import json
import sys

from filter.config_filter import config_filter
from merge.preprocessing import load_json, get_merges, get_mol, get_protein


class FilterPipeline():

    def __init__(self, merge: str, smis: list, synthons: list, fragmentA: str, fragmentB: str, proteinA: str,
                 proteinB: str, steps: list):
        self.merge = merge
        self.smis = smis
        self.synthons = synthons
        self.fragmentA = fragmentA  # filepath
        self.fragmentB = fragmentB  # filepath
        self.proteinA = proteinA  # filepath
        self.proteinB = proteinB  # filepath
        self.steps = steps  # record the filter steps to use in the pipeline
        self.mols = None
        self.results = None  # record the results for each filter (True for pass; False for fail)
        self.failures = {}  # record the SMILES that failed and at which filter

    def _remove_failed(self, filter_name: str):
        """
        Remove the failed molecules from the list of SMILES/synthons/mols to pass into the next filter.
        Record the failed SMILES and at which filter they failed in a dictionary.

        :param filter_name: name of the filter just run
        :type filter_name: str
        """
        for i, res in reversed(list(enumerate(self.results))):
            if not res:  # if merge failed the filter, then remove from data
                smi = self.smis[i]
                self.smis.pop(i)
                self.synthons.pop(i)
                if self.mols:  # if mols have been generated, remove failed merge
                    self.mols.pop(i)

                # record the failed smile in self.failures dict
                self.failures[smi] = filter_name

    def execute_pipeline(self):
        """
        For each step in the pipeline, the function retrieves the module/class and runs the filter. Each filter
        returns a list of results (True or False) and a list of molecules (returns if None if the filter does not
        generate any conformations) - the list of molecules is then passed into the next filter. Prints how many
        molecules were removed by each filter.
        """
        for step in self.steps:  # step specifies the class name for the filter
            n_smi = len(self.smis)
            module = config_filter.PIPELINE_DICT[step]  # retrieve module name
            cls = getattr(importlib.import_module(module), step)  # retrieve the class
            filter = cls(self.smis, self.synthons, self.fragmentA, self.fragmentB, self.proteinA, self.proteinB,
                         self.merge, self.mols)  # initialise the filter
            self.results, self.mols = filter.filter_all()  # run the filter
            self._remove_failed(step)  # remove failed mols from the smiles/synthon/mol lists

            n_removed = n_smi - len(self.smis)
            print(f'{n_removed} mols removed by {step}. {len(self.smis)} mols remaining.')

    def return_results(self):
        """
        Return the results after filtering.

        :return: list of filtered SMILES; dictionary with SMILES as key and failed filter as the value
        :rtype: Tuple[list, dict]
        """
        return self.smis, self.failures


def parse_args(args):
    # get file containing merges and all fragment/protein files
    parser = argparse.ArgumentParser()
    parser.add_argument('-f', '--merge_file', help='the json file containing the merges')
    parser.add_argument('-m', '--merge', help='the name of the merge, e.g. x0107_0A_x0434_0A')
    parser.add_argument('-t', '--target', help='the name of the target, e.g. nsp13 or Mpro')
    parser.add_argument('-o', '--output_dir', help='the directory to write the filtered files to')
    return parser.parse_args(args)


def main():
    args = parse_args(sys.argv[1:])
    merge = args.merge

    # open json file containing merges
    merges_dict = load_json(args.merge_file)
    synthons, smiles = get_merges(merges_dict)
    print("Number of smiles: %d"%len(smiles))

    # load fragments and proteins
    fragments = args.merge.split('-')
    fA = fragments[0]
    fB = fragments[1]
    fragmentA = get_mol(args.target, fA, config_filter.FRAGALYSIS_DATA_DIR)
    fragmentB = get_mol(args.target, fB, config_filter.FRAGALYSIS_DATA_DIR)
    proteinA = get_protein(args.target, fA, config_filter.FRAGALYSIS_DATA_DIR)
    proteinB = get_protein(args.target, fB, config_filter.FRAGALYSIS_DATA_DIR)
    steps = config_filter.PIPELINE

    # execute the pipeline
    pipeline = FilterPipeline(merge, smiles, synthons, fragmentA, fragmentB, proteinA, proteinB, steps)
    pipeline.execute_pipeline()
    print(f'{len(pipeline.smis)} mols after filtering.')
    results, failures = pipeline.return_results()

    # save in json files
    filtered_fname = merge + '_filtered.json'
    filtered_fpath = os.path.join(args.output_dir, filtered_fname)

    failed_fname = merge + '_failures.json'
    failed_fpath = os.path.join(args.output_dir, failed_fname)

    with open(filtered_fpath, 'w') as f:
        json.dump(results, f)

    with open(failed_fpath, 'w') as f:
        json.dump(failures, f)


if __name__ == "__main__":
    main()
