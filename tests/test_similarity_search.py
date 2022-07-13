import os
import unittest

import numpy as np
from rdkit import DataStructs

class TestOverlapFilter(unittest.TestCase):
    """Tests the overlap filter functions"""

    def test_morgan_fp(self):
        """Tests the function calculate distances correctly"""
        smi = "CCCOC"
        from similaritySearch.compute_fingerprints import get_fingerprint
        import similaritySearch.similaritySearchConfig as config
        config.FINGERPRINT_TYPE="morgan"
        fp = get_fingerprint(smi)
        print(fp)

    def test_ph4_fp(self):
        """Tests the function calculate distances correctly"""
        smi = "c1ccccc1C(=O)C(=O)c2ccccc2"
        from similaritySearch.compute_fingerprints import get_fingerprint
        import similaritySearch.similaritySearchConfig as config
        config.FINGERPRINT_TYPE="pharmacophore"
        fp = get_fingerprint(smi)
        print(fp)

    def test_fp_fraction_of_query_on_bits(self):
        import json
        with open(os.path.join(os.path.dirname(__file__),
                               "test_data/for_similarity/bits_morgan_2024_radidus_2.json")) as f:
            fps = json.load(f)
        print(fps)

        from similaritySearch.similarity_searcher_search_onePartition import search_fps_matrix
        from similaritySearch.compute_fingerprints import get_fingerPrint_as_npBool
        import similaritySearch.similaritySearchConfig as config

        config.FINGERPRINT_TYPE = "morgan"
        # database_dir = os.path.expanduser("~/oxford/enamine/fingerprints_db")
        database_dir = os.path.expanduser("~/localData/catalogue_search/enamine/fingerprints_db/")
        query_fps = np.zeros((3,2048), dtype=np.bool)
        dataset_names = ["CBB_all", "CBB_lateral", "CBB_diving"]
        for i, name in enumerate(dataset_names):
            query_fps[i, fps[name]] = 1

        smi = "COCCNC(=O)N1CCN(C(=O)c2cccs2)CC1C" #"COCCNC(=O)N1CCN(C(=O)c2ccco2)CC1" #"CC(C)(C)OC(=O)NC=1C=CC=C(CNC(=O)N2CCN(CC2)C(=O)C3=CC=CO3)C1"
        from similaritySearch.compute_metrics import fraction_of_query_on_bits
        target_fps = get_fingerPrint_as_npBool(smi)
        for i, name in enumerate(dataset_names):
            frac = fraction_of_query_on_bits(query_fps[i,:], target_fps)
            print(f"{name} {frac}")

    def _test_create_db(self): # _ to skip this test
        import similaritySearch.similaritySearchConfig as config
        config.FINGERPRINT_TYPE = "pharmacophore" #"morgan"
        from similaritySearch.create_db import create_db_from_multiple_files
        cxsmiles_dir = os.path.expanduser("~/oxford/enamine/cxsmiles")
        fp_outdir = os.path.expanduser("~/oxford/enamine/fingerprints_db")
        create_db_from_multiple_files(cxsmiles_dir, fp_outdir)