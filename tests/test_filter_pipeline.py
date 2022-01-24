"""Test the filter pipeline script"""
import argparse
import os
import unittest
import importlib
import json

from filter.filter_pipeline import parse_args, FilterPipeline
from merge.preprocessing import load_json, get_merges, get_mol, get_protein
from filter.config_filter import config_filter

filtered_path = 'tests/test_output/x0034_0B-x0311_0B_filtered.json'
failed_path = 'tests/test_output/x0034_0B-x0311_0B_failures.json'

class TestFilterPipeline(unittest.TestCase):
    """Tests the filter pipeline is working"""

    def test_parse_args(self):
        """
        Tests that the arg parse is working
        """
        parser = parse_args(['-f', 'tests/test_data/x0034_0B_x0311_0B.json',
                             '-m', 'x0034_0B-x0311_0B',
                             '-t', 'nsp13',
                             '-o', 'tests/test_data'])
        self.assertEqual(parser.merge_file, 'tests/test_data/x0034_0B_x0311_0B.json')
        self.assertEqual(parser.merge, 'x0034_0B-x0311_0B')
        self.assertEqual(parser.target, 'nsp13')
        self.assertEqual(parser.output_dir, 'tests/test_data')

    def test_filtering(self):
        """
        Tests a basic filter pipeline and the results - smiles undergo the descriptor, embedding and overlap filters.
        """
        json_file = 'tests/test_data/x0034_0B_x0311_0B.json'
        merge = 'x0034_0B-x0311_0B'
        target = 'nsp13'
        output_dir = 'tests/test_output'

        # open json file containing merges
        merges_dict = load_json(json_file)
        synthons, smiles = get_merges(merges_dict)
        print("Number of smiles: %d" % len(smiles))

        # load fragments and proteins
        fragments = merge.split('-')
        fA = fragments[0]
        fB = fragments[1]
        fragmentA = get_mol(target, fA, config_filter.FRAGALYSIS_DATA_DIR)
        fragmentB = get_mol(target, fB, config_filter.FRAGALYSIS_DATA_DIR)
        proteinA = get_protein(target, fA, config_filter.FRAGALYSIS_DATA_DIR)
        proteinB = get_protein(target, fB, config_filter.FRAGALYSIS_DATA_DIR)
        filter_steps = ['DescriptorFilter', 'EmbeddingFilter', 'OverlapFilter']
        score_steps = []

        # execute the pipeline
        pipeline = FilterPipeline(merge, smiles, synthons, fragmentA, fragmentB, proteinA, proteinB, filter_steps,
                                  score_steps)
        pipeline.execute_pipeline()
        results, failures = pipeline.return_results()

        # check if dictionary (exact results will change depending on parameters in config file)
        self.assertIsInstance(results, dict)
        self.assertIsInstance(failures, dict)


    def test_scoring(self):
        """
        Tests a basic filter pipeline and the results - smiles undergo the descriptor, embedding and overlap filters.
        Fragmenstein will fail - need to find a molecule that will pass...
        """
        merge = 'x107-0A-x0678-0A'
        smiles = ['NC(=O)CN1CCC2(C1)CC1(C2)OCCO1']
        synthons = ['NC(=O)C[Xe]']
        frag_dir = os.path.join('tests', 'test_Fragalysis')
        fragmentA = get_mol('Mpro', 'x0107_0A', frag_dir)
        fragmentB = get_mol('Mpro', 'x0678_0A', frag_dir)
        proteinA = get_protein('Mpro', 'x0107_0A', frag_dir)
        proteinB = get_protein('Mpro', 'x0678_0A', frag_dir)
        filter_steps = ['FragmensteinFilter']
        score_steps = ['IfpScore']
        output_dir = 'tests/test_output'

        # execute the pipeline
        pipeline = FilterPipeline(merge, smiles, synthons, fragmentA, fragmentB, proteinA, proteinB, filter_steps,
                                  score_steps)
        pipeline.execute_pipeline()
        results, failures = pipeline.return_results()

        print(results)
        print(failures)

        self.assertIsInstance(results, dict)
        self.assertIsInstance(failures, dict)


if __name__ == '__main__':
    unittest.main()
