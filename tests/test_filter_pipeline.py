"""Test the filter pipeline script"""

import json
import os
import shutil
import unittest

from filter.config_filter import config_filter
from filter.filter_pipeline import (FilterPipeline, create_directories,
                                    parse_args)
from utils.utils import get_merges, get_mol, get_protein, load_json

filtered_path = "tests/test_output/x0034_0B-x0311_0B_filtered.json"
failed_path = "tests/test_output/x0034_0B-x0311_0B_failures.json"


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

    def test_create_directories(self):
        working_dir = os.path.join("tests", "test_working")
        output_dir = os.path.join("tests", "test_output")
        pair = "x1234_0B_x4567_0B"
        target = "target"
        create_directories(target, pair, working_dir, output_dir)
        o_res = os.path.isdir(os.path.join(output_dir, target, pair))
        self.assertTrue(o_res)
        remove_files()

    def test_check_run(self):
        pipeline = FilterPipeline(
            "x1234_0B_x4567_0B",
            ["smiles1", "smiles2", "smiles3"],
            ["synthon1", "synthon2", "synthon3"],
            "fragmentA_fpath",
            "fragmentB_fpath",
            "proteinA_fpath",
            "proteinB_fpath",
            ["filter1", "filter2"],
            ["score1"],
            "target",
            os.path.join("tests", "test_working"),
            os.path.join("tests", "test_output"),
        )
        create_directories(
            "target",
            "x1234_0B_x4567_0B",
            os.path.join("tests", "test_working"),
            os.path.join("tests", "test_output"),
        )
        # create dummy file
        json_dict = {
            "x1234_0B_x4567_0B_1": {
                "pair": "x1234_0B_x4567_0B",
                "smiles": "smiles2",
                "synthon": "synthon2",
                "failed_filter": "filter1",
            }
        }
        print("£££")
        failed_fpath = os.path.join(
            pipeline.pair_output_dir, f"{pipeline.merge}_failures.json"
        )
        print(failed_fpath)
        with open(failed_fpath, "w") as f:
            json.dump(json_dict, f)

        pipeline.check_run()
        smis = pipeline.smis

        self.assertEqual(smis, ["smiles1", "smiles3"])
        remove_files()

    def test_parse_args(self):
        """
        Tests that the arg parse is working
        """
        parser = parse_args(
            [
                "-f",
                os.path.join("tests", "test_data", "x0034_0B_x0311_0B.json"),
                "-a",
                "x0034_0B",
                "-b",
                "x0311_0B",
                "-t",
                "nsp13",
                "-o",
                os.path.join("tests", "test_output"),
                "-w",
                os.path.join("tests", "test_working"),
            ]
        )
        self.assertEqual(
            parser.merge_file,
            os.path.join("tests", "test_data", "x0034_0B_x0311_0B.json"),
        )
        self.assertEqual(parser.fragmentA, "x0034_0B")
        self.assertEqual(parser.fragmentB, "x0311_0B")
        self.assertEqual(parser.target, "nsp13")
        self.assertEqual(parser.output_dir, os.path.join("tests", "test_output"))
        self.assertEqual(parser.working_dir, os.path.join("tests", "test_working"))

    def test_filtering_1(self):
        """
        Tests a basic filter pipeline and the results - smiles undergo the descriptor, embedding and overlap filters.
        """
        json_file = os.path.join("tests", "test_data", "x0034_0B_x0311_0B.json")
        target = "nsp13"
        fA = "x0034_0B"
        fB = "x0311_0B"
        merge = fA + "_" + fB
        merge = merge.replace("_", "-")
        output_dir = os.path.join("tests", "test_output")
        working_dir = os.path.join("tests", "test_working")
        create_directories("nsp13", "x0034-0B-x0311-0B", working_dir, output_dir)

        # open json file containing merges
        merges_dict = load_json(json_file)
        synthons, smiles = get_merges(merges_dict)
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
            "ExpansionFilter",
            "EmbeddingFilter",
            "OverlapFilter",
        ]
        score_steps = []

        # execute the pipeline
        pipeline = FilterPipeline(
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
            working_dir,
            output_dir,
        )
        pipeline.execute_pipeline()
        results, failures = pipeline.return_results()
        # check if dictionary (exact results will change depending on parameters in config file)
        self.assertIsInstance(results, dict)
        self.assertIsInstance(failures, dict)
        remove_files()

    def test_filtering_2(self):
        """
        Tests a basic filter pipeline and the results - smiles undergo the descriptor, embedding and overlap filters.
        Fragmenstein will fail - need to find a molecule that will pass...
        """
        merge = "x0107-0A-x0678-0A"
        smiles = ["NC(=O)CN1CCC2(C1)CC1(C2)OCCO1"]
        synthons = ["NC(=O)C[Xe]"]
        frag_dir = os.path.join("tests", "test_Fragalysis")
        fragmentA = get_mol("Mpro", "x0107_0A", False, frag_dir)
        fragmentB = get_mol("Mpro", "x0678_0A", False, frag_dir)
        proteinA = get_protein("Mpro", "x0107_0A", False, frag_dir)
        proteinB = get_protein("Mpro", "x0678_0A", False, frag_dir)
        target = "Mpro"
        filter_steps = ["FragmensteinFilter"]
        score_steps = ["IfpScore"]
        output_dir = "tests/test_output"
        working_dir = "tests/test_working"

        create_directories("Mpro", "x0107-0A-x0678-0A", working_dir, output_dir)
        # execute the pipeline
        pipeline = FilterPipeline(
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
            working_dir,
            output_dir,
        )
        pipeline.execute_pipeline()
        results, failures = pipeline.return_results()

        self.assertIsInstance(results, dict)
        self.assertIsInstance(failures, dict)
        remove_files()

    def test_scoring(self):
        """
        Tests a basic filter pipeline and the results - smiles undergo the descriptor, embedding and overlap filters.
        Fragmenstein will fail - need to find a molecule that will pass...
        """
        merge = "x0991-0A-x0072-0A"
        data = load_json(os.path.join("tests", "test_data", "x0991_0A_x0072_0A.json"))
        smiles = data["smiles"][6:8]
        synthons = data["synthons"][6:8]
        fragmentA = get_mol("Mpro", "x0991_0A", False)
        fragmentB = get_mol("Mpro", "x0072_0A", False)
        proteinA = get_protein("Mpro", "x0991_0A", False)
        proteinB = get_protein("Mpro", "x0072_0A", False)
        target = "Mpro"
        filter_steps = ["FragmensteinFilter"]
        score_steps = ["IfpScore"]
        output_dir = "tests/test_output"
        working_dir = "tests/test_working"

        create_directories("Mpro", "x0991-0A-x0072-0A", working_dir, output_dir)
        # execute the pipeline
        pipeline = FilterPipeline(
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
            working_dir,
            output_dir,
        )
        pipeline.execute_pipeline()
        results, failures = pipeline.return_results()

        self.assertIsInstance(results, dict)
        self.assertIsInstance(failures, dict)

        remove_files()

    def test_updated_pipeline(self):
        """
        Tests a basic filter pipeline and the results - smiles undergo the descriptor, embedding and overlap filters.
        Fragmenstein will fail - need to find a molecule that will pass...
        """
        merge = "x0991-0A-x0072-0A"
        data = load_json(os.path.join("tests", "test_data", "x0991_0A_x0072_0A.json"))
        smiles = data["smiles"][6:8]
        synthons = data["synthons"][6:8]
        fragmentA = get_mol("Mpro", "x0991_0A", False)
        fragmentB = get_mol("Mpro", "x0072_0A", False)
        proteinA = get_protein("Mpro", "x0991_0A", False)
        proteinB = get_protein("Mpro", "x0072_0A", False)
        target = "Mpro"
        filter_steps = ["FragmensteinFilter", "EnergyFilter"]
        score_steps = ["PlipIfpScore", "SuCOSScore", "ElaboratabilityScore"]
        output_dir = "tests/test_output"
        working_dir = "tests/test_working"

        merge_dir = create_directories("Mpro", "x0991-0A-x0072-0A", working_dir, output_dir)
        # execute the pipeline
        pipeline = FilterPipeline(
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
