import os

from similaritySearch.similarity_searcher_search_onePartition import search_smi_list

DB_DIR = os.path.expanduser("/tmp/fingerprints_db_debug")

if __name__ == "__main__":
    query_smi_list = [["CCN(CC1=CC=C(F)C(F)=C1)CC(O)COCCCO"], ["CCCOCCCO"], ["CCCCCCCNCOC"]]
    results = search_smi_list(query_smi_list, DB_DIR, n_hits_per_smi=10,
                              output_name=None, backend="numba", metric="Fp_bits_frequency2D", n_cpus=1, verbose=True)
    print(results)