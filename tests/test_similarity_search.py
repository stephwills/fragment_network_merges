import json
import os
import tempfile
import unittest

import numpy as np

from fragment_network_merges.similaritySearch import search_smi_list


class TestSimilaritySearch(unittest.TestCase):
    """Tests the overlap filter functions"""

    DB_DIR = os.path.expanduser("~/oxford/enamine/fingerprints_db")

    def test_morgan_fp(self):
        """Tests the function calculate distances correctly"""
        smi = "CCCOC"
        from fragment_network_merges.similaritySearch import get_fingerprint
        import fragment_network_merges.similaritySearch.similaritySearchConfig as config
        config.FINGERPRINT_TYPE="morgan"
        fp = get_fingerprint(smi)
        print(fp)

    def test_ph4_fp(self):
        """Tests the function calculate distances correctly"""
        smi = "c1ccccc1C(=O)C(=O)c2ccccc2"
        from fragment_network_merges.similaritySearch import get_fingerprint
        import fragment_network_merges.similaritySearch.similaritySearchConfig as config
        config.FINGERPRINT_TYPE="pharmacophore"
        fp = get_fingerprint(smi)
        print(fp)

    def test_fp_fraction_of_query_on_bits(self):
        import json
        with open(os.path.join(os.path.dirname(__file__),
                               "test_data/for_similarity/bits_morgan_2024_radidus_2.json")) as f:
            fps = json.load(f)
        print(fps)

        from fragment_network_merges.similaritySearch import get_fingerPrint_as_npBool
        import fragment_network_merges.similaritySearch.similaritySearchConfig as config

        config.FINGERPRINT_TYPE = "morgan"
        # database_dir = os.path.expanduser("~/oxford/enamine/fingerprints_db")
        database_dir = os.path.expanduser("~/localData/catalogue_search/enamine/fingerprints_db/")
        query_fps = np.zeros((3,2048), dtype=np.bool)
        dataset_names = ["CBB_all", "CBB_lateral", "CBB_diving"]
        for i, name in enumerate(dataset_names):
            query_fps[i, fps[name]] = 1

        smi = "COCCNC(=O)N1CCN(C(=O)c2cccs2)CC1C" #"COCCNC(=O)N1CCN(C(=O)c2ccco2)CC1" #"CC(C)(C)OC(=O)NC=1C=CC=C(CNC(=O)N2CCN(CC2)C(=O)C3=CC=CO3)C1"
        from fragment_network_merges.similaritySearch.compute_metrics import fraction_of_query_on_bits
        target_fps = get_fingerPrint_as_npBool(smi)
        for i, name in enumerate(dataset_names):
            frac = fraction_of_query_on_bits(query_fps[i,:], target_fps)
            print(f"{name} {frac}")

    def test_search_one_partition(self):
        query_smi_list= [["CCCOCCCO"], ["CCCCCCCNCOC"]]
        results = search_smi_list(query_smi_list, self.DB_DIR, n_hits_per_smi=10,
                        output_name=None, backend="numba", metric="Tanimoto", n_cpus=1, verbose=True)
        print(results)

    def test_search_one_partition_tanimotoW(self):
        query_smi_list= [["C=CC(C)C(NC(=O)CC1=CC=C(S(N)(=O)=O)S1)C(=O)OC"], ["CCCCCCCNCOC"]]
        results = search_smi_list(query_smi_list, self.DB_DIR, n_hits_per_smi=10,
                        output_name=None, backend="numba", metric="TanimotoW", n_cpus=1, verbose=True)
        print(results)

    def test_search_one_partition_Fp_bits_frequency(self):
        query_smi_list= [["CCN(CC1=CC=C(F)C(F)=C1)CC(O)COCCCO"], ["CCCOCCCO"], ["CCCCCCCNCOC"]]
        results = search_smi_list(query_smi_list, self.DB_DIR, n_hits_per_smi=10,
                        output_name=None, backend="numba", metric="Fp_bits_frequency", n_cpus=1, verbose=True)
        print(results)
        "(\d+\.\d+\s+\w+)\s+(.+)"

    def test_search_one_partition_Fp_bits_frequency2D(self):
        query_smi_list= [["CCN(CC1=CC=C(F)C(F)=C1)CC(O)COCCCO"], ["CCCOCCCO"], ["CCCCCCCNCOC"]]
        results = search_smi_list(query_smi_list, self.DB_DIR, n_hits_per_smi=10,
                        output_name=None, backend="numba", metric="Fp_bits_frequency2D", n_cpus=1, verbose=True)
        print(results)

    def test_search_one_partition_Fp_bits_frequency_combiningFps(self):
        query_smi_list= [["CCN(CC1=CC=C(F)C(F)=C1)CC(O)COCCCO", "FCCCCCCCCCCCNCOC", "CN1C=NC2=C1C(=O)N(C(=O)N2C)C",
                          "CN1CCOC(C)(C(=O)NCc2ccccn2)C1", "S=C(Nc1ccccc1)N1CCN(Cc2ccccn2)CC1"]]
        results = search_smi_list(query_smi_list, self.DB_DIR, n_hits_per_smi=10, conserved_common_bits_fraction=0.5,
                        output_name=None, backend="numba", metric="Fp_bits_frequency", n_cpus=1, verbose=True)
        print(results)


    def test_numba_logexpsum(self):
        in_ = np.random.rand(10, 8)
        from fragment_network_merges.similaritySearch.compute_metrics import numba_logsumexp_stable
        out = numba_logsumexp_stable(in_)
        from scipy.special import logsumexp
        out2 = logsumexp(in_, axis=1, return_sign=False)
        self.assertTrue(np.isclose(out, out2).all())

    def test_cli_one_partion(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            cmd =f'python ' \
                 f'-m similaritySearch.similarity_searcher_search_onePartition -d {self.DB_DIR} ' \
                 f'--verbose --n_cpus 2 -o {tmpdir}/first_partition.json -'
            import subprocess
            print(cmd)
            subprocess.run(cmd, shell=True, input="COC(C)CS(=O)(=O)NCC1=C(C)C(C)=C(C)C(C)=C1C\nCOC\nC1C(C(C(C(O1)O)O)O)O".encode())
            with open(f'{tmpdir}/first_partition.json') as f:
                data = json.load(f)
                self.assertAlmostEqual(data["COC(C)CS(=O)(=O)NCC1=C(C)C(C)=C(C)C(C)=C1C"][0][0], 1.0)

    def _test_create_db(self): # _ to skip this test
        import fragment_network_merges.similaritySearch.similaritySearchConfig as config
        config.FINGERPRINT_TYPE = "pharmacophore" #"morgan"
        from fragment_network_merges.similaritySearch import create_db_from_multiple_files
        cxsmiles_dir = os.path.expanduser("~/oxford/enamine/cxsmiles")
        fp_outdir = self.DB_DIR
        create_db_from_multiple_files(cxsmiles_dir, fp_outdir)