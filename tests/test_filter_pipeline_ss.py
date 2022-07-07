"""Test the filter pipeline script"""

import json
import os
import shutil
import unittest

from filter.config_filter import config_filter
from filter.filter_pipeline_ss import FilterPipelineSS, create_directories, parse_args
from utils.utils import get_mol, get_protein, load_json, get_merges


def remove_files():
    w_dir = os.path.join("tests", "test_working")
    o_dir = os.path.join("tests", "test_output")
    for dir in [w_dir, o_dir]:
        for file in os.listdir(dir):
            path = os.path.join(dir, file)
            try:
                shutil.rmtree(path)
            except OSError:
                os.remove(path)


class TestFilterPipeline(unittest.TestCase):
    """Tests the filter pipeline is working"""

    def test_filtering_1(self):
        """
        Tests a basic filter pipeline and the results - smiles undergo the descriptor, embedding and overlap filters.
        """
        json_file = os.path.join("tests", "test_data", "x0034_0B_x0176_0B_ss.json")
        target = "nsp13"
        fA = "x0034_0B"
        fB = "x0176_0B"
        merge = fA + "_" + fB
        merge = merge.replace('_', '-')
        output_dir = os.path.join("tests", "test_output")
        working_dir = os.path.join("tests", "test_working")
        merge_dir = create_directories("nsp13", "x0034-0B-x0176-0B", working_dir, output_dir)

        # open json file containing merges
        smiles = load_json(json_file)
        synthons = None
        print("Number of smiles: %d" % len(smiles))

        # load fragments and proteins
        fragmentA = get_mol(
            target, fA, fragalysis_dir=config_filter.FRAGALYSIS_DATA_DIR
        )
        fragmentB = get_mol(
            target, fB, fragalysis_dir=config_filter.FRAGALYSIS_DATA_DIR
        )
        proteinA = get_protein(
            target, fA, fragalysis_dir=config_filter.FRAGALYSIS_DATA_DIR
        )
        proteinB = get_protein(
            target, fB, fragalysis_dir=config_filter.FRAGALYSIS_DATA_DIR
        )
        filter_steps = ["DescriptorFilter", "EmbeddingFilterSS", "OverlapFilter"]
        score_steps = []

        # execute the pipeline
        pipeline = FilterPipelineSS(
            merge,
            smiles,
            synthons,
            fragmentA,
            fragmentB,
            proteinA,
            proteinB,
            filter_steps,
            score_steps,
            target,
            merge_dir,
            working_dir,
            output_dir,
        )
        pipeline.execute_pipeline()
        results, failures = pipeline.return_results()
        # check if dictionary (exact results will change depending on parameters in config file)
        self.assertIsInstance(results, dict)
        self.assertIsInstance(failures, dict)
        remove_files()

    def test_scoring(self):
        """
        Tests a basic filter pipeline and the results - smiles undergo the descriptor, embedding and overlap filters.
        Fragmenstein will fail - need to find a molecule that will pass...
        """
        json_file = os.path.join("tests", "test_data", "x0034_0B_x0176_0B_ss.json")
        target = "nsp13"
        fA = "x0034_0B"
        fB = "x0176_0B"
        merge = fA + "_" + fB
        merge = merge.replace("_", "-")
        output_dir = os.path.join("tests", "test_output")
        working_dir = os.path.join("tests", "test_working")
        merge_dir = create_directories("nsp13", "x0034-0B-x0176-0B", working_dir, output_dir)

        # open json file containing merges
        smiles = load_json(json_file)
        synthons = None
        print("Number of smiles: %d" % len(smiles))

        # load fragments and proteins
        fragmentA = get_mol(
            target, fA, fragalysis_dir=config_filter.FRAGALYSIS_DATA_DIR
        )
        fragmentB = get_mol(
            target, fB, fragalysis_dir=config_filter.FRAGALYSIS_DATA_DIR
        )
        proteinA = get_protein(
            target, fA, fragalysis_dir=config_filter.FRAGALYSIS_DATA_DIR
        )
        proteinB = get_protein(
            target, fB, fragalysis_dir=config_filter.FRAGALYSIS_DATA_DIR
        )

        filter_steps = [
            "DescriptorFilter",
            "EmbeddingFilterSS",
            "OverlapFilter",
            "FragmensteinFilter",
        ]
        score_steps = ["IfpScore"]

        # execute the pipeline
        pipeline = FilterPipelineSS(
            merge,
            smiles,
            synthons,
            fragmentA,
            fragmentB,
            proteinA,
            proteinB,
            filter_steps,
            score_steps,
            target,
            merge_dir,
            working_dir,
            output_dir,
        )
        pipeline.execute_pipeline()
        results, failures = pipeline.return_results()

        self.assertIsInstance(results, dict)
        self.assertIsInstance(failures, dict)

        remove_files()


if __name__ == "__main__":
    unittest.main()
