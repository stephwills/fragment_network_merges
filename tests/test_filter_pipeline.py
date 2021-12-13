"""Test the filter pipeline script"""
import argparse
import os
import unittest

from filter.filter_pipeline import parse_args

filtered_path = 'tests/test_output/x0034_0B-x0311_0B_filtered.json'
failed_path = 'tests/test_output/x0034_0B-x0311_0B_failures.json'


def remove_files():
    if os.path.exists(filtered_path):
        os.remove(filtered_path)
    if os.path.exists(failed_path):
        os.remove(failed_path)
    print('Files removed')


class TestFilterPipeline(unittest.TestCase):
    """Tests the filter pipeline is working"""

    def test_parse_args(self):
        parser = parse_args(['-f', 'tests/test_data/x0034_0B_x0311_0B.json',
                             '-m', 'x0034_0B-x0311_0B',
                             '-t', 'nsp13',
                             '-o', 'tests/test_data'])
        self.assertEqual(parser.merge_file, 'tests/test_data/x0034_0B_x0311_0B.json')
        self.assertEqual(parser.merge, 'x0034_0B-x0311_0B')
        self.assertEqual(parser.target, 'nsp13')
        self.assertEqual(parser.output_dir, 'tests/test_data')

    def test_main(self):
        os.system('python filter/filter_pipeline.py -f tests/test_data/x0034_0B_x0311_0B.json -m x0034_0B-x0311_0B -t nsp13 -o tests/test_output')
        self.assertTrue(os.path.exists(filtered_path))
        self.assertTrue(os.path.exists(failed_path))
        remove_files()


if __name__ == '__main__':
    unittest.main()
