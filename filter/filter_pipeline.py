"""Run the filtering pipeline"""
import argparse
import importlib
import os
import json
import sys
import shutil
import time
import tempfile

from filter.config_filter import config_filter
from utils.utils import load_json, get_merges, get_mol, get_protein


def create_directories(
    target: str,
    pair: str,
    working_dir: str = config_filter.WORKING_DIR,
    output_dir: str = config_filter.OUTPUT_DIR,
):
    """
    Create directories to write files to. Write intermediate files to a temporary location. Sub-directories created
    for each pair.

    :param pair: name of pair
    :type pair: str
    :param working_dir: working directory path
    :type working_dir: str
    :param output_dir: output directory path
    :type output_dir: str
    """
    tmpdir = tempfile.mkdtemp(dir=working_dir, prefix='tempfiles_')
    merge_dir = os.path.join(tmpdir, target, pair)

    if not os.path.exists(os.path.join(tmpdir, target)):
        os.mkdir(os.path.join(tmpdir, target))

    if not os.path.exists(merge_dir):
        os.mkdir(merge_dir)

    if not os.path.exists(os.path.join(output_dir, target)):
        os.mkdir(os.path.join(output_dir, target))

    if not os.path.exists(os.path.join(output_dir, target, pair)):
        os.mkdir(os.path.join(output_dir, target, pair))

    return merge_dir


class FilterPipeline:
    def __init__(
        self,
        merge: str,
        smis: list,
        synthons,  # None if sim search; list if frag net
        fragmentA: str,
        fragmentB: str,
        proteinA: str,
        proteinB: str,
        filter_steps: list,
        score_steps: list,
        target: str,
        merge_dir: str,
        working_dir: str = config_filter.WORKING_DIR,
        output_dir: str = config_filter.OUTPUT_DIR,
    ):
        self.merge = merge
        self.smis = smis
        self.synthons = synthons
        self.fragmentA = fragmentA  # filepath
        self.fragmentB = fragmentB  # filepath
        self.proteinA = proteinA  # filepath
        self.proteinB = proteinB  # filepath
        self.filter_steps = (
            filter_steps  # record the filter steps to use in the pipeline
        )
        self.score_steps = (
            score_steps  # record the scoring metrics to score the filtered merges
        )
        self.mols = None
        self.results = (
            None  # record the results for each filter (True for pass; False for fail)
        )
        self.failures = {}  # record the SMILES that failed and at which filter
        self.timings = {}
        self.score_dict = {}
        self.target = target

        # file saving
        self.pair_working_dir = merge_dir  # the pair dir within the working directory
        self.pair_output_dir = os.path.join(output_dir, target, merge)  # the pair dir within the output directory
        self.failed_fpath = os.path.join(
            self.pair_working_dir, f"{self.merge}_failures.json"
        )
        self.timings_fpath = os.path.join(
            self.pair_working_dir, f"{self.merge}_timings.json"
        )

        # get unique merge name - will be used for naming output files
        self.names = [self.merge + "_" + str(i) for i, _ in enumerate(self.smis)]

        # store filepaths generated by the filtering pipeline
        self.mol_files = None
        self.holo_files = None
        self.apo_files = None

    def check_run(self):
        """
        Check if any of the SMILES have already been run - this happens.
        """
        output_failed_fpath = os.path.join(
            self.pair_output_dir, f"{self.merge}_failures.json"
        )
        print('$$$')
        print(output_failed_fpath)
        if os.path.exists(output_failed_fpath):
            self.failures = load_json(output_failed_fpath)
            failed_merges = list(self.failures.keys())

            for i, name in reversed(list(enumerate(self.names))):
                if name in failed_merges:
                    self.smis.pop(i)
                    self.names.pop(i)
                    if self.synthons:
                        self.synthons.pop(i)

            print(
                f"{len(failed_merges)} merges have already been filtered. {len(self.smis)} merges remaining."
            )

        if os.path.exists(self.timings_fpath):
            self.timings = load_json(self.timings_fpath)

    def _remove_failed(self, filter_name: str, working_dir=config_filter.WORKING_DIR):
        """
        Remove the failed molecules from the list of SMILES/synthons/mols to pass into the next filter.
        Record the failed SMILES and at which filter they failed in a dictionary.
        Write the failed molecules after each filter applied.

        :param filter_name: name of the filter just run
        :type filter_name: str
        """
        for i, res in reversed(list(enumerate(self.results))):
            if not res:  # if merge failed the filter, then remove from data
                smi = self.smis[i]
                if self.synthons:
                    synthon = self.synthons[i]
                name = self.names[i]

                # record the failed smile in self.failures dict
                if self.synthons:
                    inner_dict = {
                        "pair": self.merge,
                        "smiles": smi,
                        "synthon": synthon,
                        "failed_filter": filter_name,
                    }
                else:
                    inner_dict = {
                        "pair": self.merge,
                        "smiles": smi,
                        "failed_filter": filter_name,
                    }
                self.failures[name] = inner_dict
                self.smis.pop(i)
                if self.synthons:
                    self.synthons.pop(i)
                self.names.pop(i)
                if self.mols:  # if mols have been generated, remove failed merge
                    self.mols.pop(i)
                if self.mol_files:
                    self.mol_files.pop(i)
                if self.apo_files:
                    self.apo_files.pop(i)
                if self.holo_files:
                    self.holo_files.pop(i)

    def _write_log(self, failed_fpath=None, timings_fpath=None):
        """
        Write failure dict to file (help prevent re-running SMILES later)
        """
        if os.path.exists(self.pair_working_dir):
            with open(self.failed_fpath, "w") as f:
                json.dump(self.failures, f)

            with open(self.timings_fpath, "w") as f:
                json.dump(self.timings, f)

        else:
            if failed_fpath:
                with open(failed_fpath, "w") as f:
                    json.dump(self.failures, f)

            if timings_fpath:
                with open(timings_fpath, "w") as f:
                    json.dump(self.timings, f)

    def _move_output(self):
        """
        After each filter, move the output files to output directory.
        """
        files = os.listdir(self.pair_working_dir)
        files = [os.path.join(self.pair_working_dir, f) for f in files]
        for f in files:
            if os.path.isfile(f):
                new_fpath = os.path.join(self.pair_output_dir, os.path.basename(f))
                shutil.copy(f, new_fpath)

    def execute_pipeline(self):
        """
        For each step in the pipeline, the function retrieves the module/class and runs the filter. Each filter
        returns a list of results (True or False) and a list of molecules (returns if None if the filter does not
        generate any conformations) - the list of molecules is then passed into the next filter. Prints how many
        molecules were removed by each filter.
        """
        for step in self.filter_steps:  # step specifies the class name for the filter
            start = time.time()
            n_smi = len(self.smis)
            module = config_filter.PIPELINE_DICT[step]  # retrieve module name
            cls = getattr(importlib.import_module(module), step)  # retrieve the class
            filter = cls(
                self.smis,
                self.synthons,
                self.fragmentA,
                self.fragmentB,
                self.proteinA,
                self.proteinB,
                self.merge,
                self.mols,
                self.names,
                self.pair_working_dir,
                self.pair_output_dir,
            )  # initialise the filter
            self.results, self.mols = filter.filter_all()  # run the filter
            # if the filter creates placed mol/pdb files, retrieve them (will stay as None if not)
            if filter.get_placed_files()[0]:
                self.mol_files, self.holo_files, self.apo_files = filter.get_placed_files()
            end = time.time()
            self._remove_failed(
                step
            )  # remove failed mols from the smiles/synthon/mol lists
            n_removed = n_smi - len(self.smis)
            if step not in self.timings.keys():
                timings_dict = {
                    "time": round(end - start, 2),
                    "n_run": n_smi,
                    "n_removed": n_removed,
                }  # timings to file
                self.timings[step] = timings_dict
            self._write_log()
            print(
                f"{n_removed} mols removed by {step}. {len(self.smis)} mols remaining."
            )
            self._move_output()

        if (
            len(self.smis) > 0
            and len(self.score_steps) > 0
            and self.mol_files
            or self.holo_files
            or self.apo_files
        ):
            start = time.time()
            print("Beginning scoring")

            for step in self.score_steps:
                module = config_filter.PIPELINE_DICT[step]  # retrieve module name
                cls = getattr(
                    importlib.import_module(module), step
                )  # retrieve the class
                score = cls(
                    self.smis,
                    self.synthons,
                    self.fragmentA,
                    self.fragmentB,
                    self.proteinA,
                    self.proteinB,
                    self.merge,
                    self.mols,
                    self.names,
                    self.pair_working_dir,
                    self.pair_output_dir,
                    self.mol_files,
                    self.holo_files,
                    self.apo_files,
                )
                self.score_dict[step] = score.score_all()
                print(f"Scoring metric {step} run on filtered merges.")
                if step not in self.timings.keys():
                    end = time.time()
                    timings_dict = {
                        "time": round(end - start, 2),
                        "n_scored": len(self.smis),
                    }  # timings to file
                    self.timings[step] = timings_dict
                self._write_log()
                self._move_output()

    def return_results(self):
        """
        Return the results after filtering.

        :return: list of filtered SMILES; dictionary with SMILES as key and failed filter as the value
        :rtype: Tuple[list, dict]
        """
        if len(self.score_dict) > 0:
            results_dict = {}
            if self.synthons:
                for i, (name, smi, synthon) in enumerate(
                    zip(self.names, self.smis, self.synthons)
                ):
                    inner_dict = {"pair": self.merge, "smiles": smi, "synthon": synthon}
                    for score in self.score_dict:
                        inner_dict[score] = self.score_dict[score][i]
                    results_dict[name] = inner_dict

            else:
                for i, (name, smi) in enumerate(
                        zip(self.names, self.smis)
                ):
                    inner_dict = {"pair": self.merge, "smiles": smi}
                    for score in self.score_dict:
                        inner_dict[score] = self.score_dict[score][i]
                    results_dict[name] = inner_dict

        else:
            results_dict = {}
            if self.synthons:
                for i, (name, smi, synthon) in enumerate(
                    zip(self.names, self.smis, self.synthons)
                ):
                    inner_dict = {"pair": self.merge, "smiles": smi, "synthon": synthon}
                    results_dict[name] = inner_dict

            else:
                for i, (name, smi) in enumerate(
                        zip(self.names, self.smis)
                ):
                    inner_dict = {"pair": self.merge, "smiles": smi}
                    results_dict[name] = inner_dict

        # save in json files
        filtered_fname = self.merge + "_filtered.json"
        filtered_fpath = os.path.join(self.pair_output_dir, filtered_fname)

        with open(filtered_fpath, "w") as f:
            json.dump(results_dict, f)

        return results_dict, self.failures


def parse_args(args):
    # get file containing merges and all fragment/protein files
    parser = argparse.ArgumentParser(
        epilog="""
    Example
    python -m filter.filter_pipeline -f x0034_0B_x0176_0B.json -m x0034_0B_x0176_0B -t nsp13 -o data/example_folder"""
    )
    parser.add_argument(
        "-f", "--merge_file", help="the json file containing the merges"
    )
    parser.add_argument(
        "-a", "--fragmentA", help="the name of fragment A, e.g. x0107_0A"
    )
    parser.add_argument(
        "-b", "--fragmentB", help="the name of fragment B, e.g. x0434_0A"
    )
    parser.add_argument(
        "-t", "--target", help="the name of the target, e.g. nsp13 or Mpro"
    )
    parser.add_argument(
        "-o", "--output_dir", help="the directory to write the filtered files to"
    )
    parser.add_argument(
        "-w", "--working_dir", help="the directory to write the intermediate files to"
    )
    parser.add_argument(
        "-s", "--sim_search", help="whether the data is similarity search data or not", action='store_true'
    )
    return parser.parse_args(args)


def main():
    args = parse_args(sys.argv[1:])
    fA = args.fragmentA
    fB = args.fragmentB
    merge = fA + "-" + fB
    merge = merge.replace("_", "-")
    merge_dir = create_directories(args.target, merge, args.working_dir, args.output_dir)

    # open json file containing merges
    if not args.sim_search:
        merges_dict = load_json(args.merge_file)
        synthons, smiles = get_merges(merges_dict)
    else:
        smiles = load_json(args.merge_file)
        synthons = None
    print("Number of smiles: %d" % len(smiles))

    # load fragments and proteins
    fragmentA = get_mol(
        args.target, fA, fragalysis_dir=config_filter.FRAGALYSIS_DATA_DIR
    )
    fragmentB = get_mol(
        args.target, fB, fragalysis_dir=config_filter.FRAGALYSIS_DATA_DIR
    )
    proteinA = get_protein(
        args.target, fA, fragalysis_dir=config_filter.FRAGALYSIS_DATA_DIR
    )
    proteinB = get_protein(
        args.target, fB, fragalysis_dir=config_filter.FRAGALYSIS_DATA_DIR
    )
    filter_steps = config_filter.FILTER_PIPELINE
    score_steps = config_filter.SCORING_PIPELINE

    # execute the pipeline
    pipeline = FilterPipeline(
        merge,
        smiles,
        synthons,  # None if sim search data
        fragmentA,
        fragmentB,
        proteinA,
        proteinB,
        filter_steps,
        score_steps,
        args.target,
        merge_dir,
        args.working_dir,
        args.output_dir,
    )
    pipeline.check_run()
    pipeline.execute_pipeline()
    print(f"{len(pipeline.smis)} mols after filtering.")
    results, _ = pipeline.return_results()

    # save in json files
    filtered_fname = merge + "_filtered.json"
    filtered_fpath = os.path.join(args.output_dir, args.target, merge, filtered_fname)

    if not os.path.exists(filtered_fpath):
        with open(filtered_fpath, "w") as f:
            json.dump(results, f)

    shutil.rmtree(merge_dir)


if __name__ == "__main__":
    main()
