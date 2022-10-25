import datetime
import io
import json
import os
import sqlite3
import time
from collections import OrderedDict

import dask.bag as db
import numba
import numpy as np
import pandas as pd

from similaritySearch.compute_fingerprints import get_fingerPrint_as_npBool, \
    decompressFingerprint_npStr
import similaritySearch.similaritySearchConfig as config
from similaritySearch.compute_metrics import jaccard_vectorized, jaccard_numba, \
    tversky_numba, fp_bits_frequency, jaccard_weighted_numba, fp_bits_frequency2D
from utils.parallelUtils import get_parallel_client

def combine_two_chunk_searches(cs1, cs2):
    '''

    :param cs1: tuple( similarities_1, ids_1)
    :param cs2: tuple( similarities_2, ids_2)
    :return:
    '''
    assert cs1[0].shape == cs2[0].shape, "Error, union for different shapes not implemented"

    sim_concat =  np.concatenate([ cs1[0], cs2[0]], axis=1 )
    idxs_concat = np.concatenate([ cs1[1], cs2[1]], axis=1 )

    to_pick_idxs = np.argsort(-sim_concat, axis = -1)[:, :cs1[0].shape[1]]
    new_sim =  -np.inf * np.ones_like(cs1[0])
    new_idxs = -1 * np.ones_like(cs1[1])
    for i in range(sim_concat.shape[0]):
        new_sim[i,:] = sim_concat[i, to_pick_idxs[i, :]]
        new_idxs[i,...] = idxs_concat[i, to_pick_idxs[i, :], :]

    # print("s1"); print(cs1[0])
    # print("s2"); print(cs2[0])
    # print("res="); print(new_sim)

    return new_sim, new_idxs

def process_chunk_using_numpy(query_fps_mat, db_many_fps_bool, n_hits_per_smi, verbose=False ):
    sim_matrix = jaccard_vectorized(query_fps_mat, db_many_fps_bool)
    max_sim_idx_per_compound = np.argsort(-sim_matrix, axis=1)[:, :n_hits_per_smi]
    n_hits_found = max_sim_idx_per_compound.shape[1]
    found_similarities = -np.inf * np.ones((query_fps_mat.shape[0], n_hits_per_smi))
    found_similarities[:, : n_hits_found] = np.take_along_axis(sim_matrix, max_sim_idx_per_compound, axis=1)

    idxs = -1 * np.ones((query_fps_mat.shape[0], n_hits_per_smi), dtype=np.int64)
    idxs[:, : n_hits_found] = max_sim_idx_per_compound
    return found_similarities, idxs, n_hits_found

@numba.jit(**config.NUMBA_kwARGS)
def process_chunk_using_numba(query_fps_mat, db_fps_mat, n_hits_per_smi, metric_fun, verbose=False):

    n_query_fps = query_fps_mat.shape[0]
    matched_similarities = np.ones((n_query_fps, n_hits_per_smi)) * -np.inf  # query_id, hit_num, similarity
    matched_ids = np.ones((n_query_fps, n_hits_per_smi), dtype=np.int64) * -1  # query_id, hit_num,  hit_id

    n_mols = db_fps_mat.shape[0]
    n_mol_processed = 0

    fp_size = db_fps_mat.shape[1]
    query_multiplicity = query_fps_mat.shape[1] // fp_size

    print_step = 10 ** int(np.ceil(np.log10((int(1e5) // n_query_fps))))
    for j in range(db_fps_mat.shape[0]):
        # if j%10==0: print(j, db_fps_mat.shape[0])
        db_binary_finPrint = db_fps_mat[j]
        n_mol_processed += 1
        if verbose and n_mol_processed % print_step == 0: print(n_mol_processed, "mols processed in block of size", n_mols)
        for i in numba.prange(query_fps_mat.shape[0]):
            query_fp = query_fps_mat[i]
            prod = 1.
            for subquery_idx in range(query_multiplicity):
                # query_sim =metric_fun(query_fp[(subquery_idx * fp_size):((subquery_idx + 1) * fp_size)], db_binary_finPrint)
                prod *= metric_fun(query_fp[(subquery_idx * fp_size):((subquery_idx + 1) * fp_size)], db_binary_finPrint)
            prod = prod ** (1. / query_multiplicity)
            # simil = prod if 0<=prod<=1 else -1
            simil = prod
            less_similar_idx = np.argmin(matched_similarities[i, :])

            if simil > matched_similarities[i, less_similar_idx]:
                matched_similarities[i, less_similar_idx] = simil
                matched_ids[i, less_similar_idx] = j

    n_elems_found = 0
    for elem in matched_similarities[0,...]:
        if elem<=-np.inf:
            break
        n_elems_found += 1
    return matched_similarities, matched_ids, n_elems_found

def process_one_subFile_by_chunks(query_fps_mat, fileNum_chunkFname, n_hits_per_smi, process_chunk_fun,
                                  n_mols_per_chunk=1000000, metric_fun= jaccard_numba, verbose=False):
    assert len(query_fps_mat.shape) > 1, "Error query_fps_mat should be of rank 2"
    file_num, chunk_fname = fileNum_chunkFname
    initial_time = datetime.datetime.now()
    if verbose: print( "[%s]: Processing %s: %s"%(initial_time.strftime("%c"), file_num, chunk_fname))
    if n_mols_per_chunk is not None:
        chunkSize = n_mols_per_chunk * (query_fps_mat.shape[1] // 8)
    else:
        raise  ValueError("n_mols_per_chunk is required")
    matched_similarities = -np.inf * np.ones((query_fps_mat.shape[0], n_hits_per_smi))
    matched_ids = -1 * np.ones(matched_similarities.shape + (2,), dtype=np.int64)
    matched_ids[..., 0] = file_num


    with open(chunk_fname, "rb") as f:
        bytes_ = f.read(chunkSize)
        i=0
        while bytes_:
            db_many_fps_bool = decompressFingerprint_npStr( bytes_ ).reshape(-1, config.FINGERPRINT_NBITS)
            found_similarities, idxs, n_hits_found = process_chunk_fun( query_fps_mat, db_many_fps_bool, n_hits_per_smi, metric_fun )
            idxs[:, : n_hits_found] += i
            idxs = np.stack([file_num*np.ones_like(idxs), idxs], axis = -1 )
            matched_similarities, matched_ids= combine_two_chunk_searches((found_similarities, idxs), (matched_similarities, matched_ids))

            i += db_many_fps_bool.shape[0]
            bytes_ = f.read(chunkSize)

    if verbose: print( "[%s]: %s: %s was computed in %s"%( datetime.datetime.now().strftime("%c"), chunk_fname, file_num, datetime.datetime.now()-initial_time))
    return matched_similarities, matched_ids

def search_smi_list(query_smi_list, database_dir, n_hits_per_smi=30, output_name=None,
                    conserved_common_bits_fraction = -1, backend="numpy",
                    metric="Tanimoto", verbose=True, n_cpus=1):
    """
    param: query_smi_list: List of lists of strings representing smiles to be queried
    """
    n_smiles_per_group = [len(elem)  for elem in query_smi_list]
    assert all([elem > 0 for elem in n_smiles_per_group]), "Error, for combine_fingerprint options, more than one" \
                                                            " smiles per group should be provided"
    if 0 < conserved_common_bits_fraction <1:
        combine_fps = True
    else:
        combine_fps = False

    starting_time = time.time()
    def compute_query_fp(qsmi_list):
        fps=[]
        for qsmi in qsmi_list:
            fp = get_fingerPrint_as_npBool(qsmi)
            if fp is None:
                return None
            else:
                fps.append(fp)
        if combine_fps:
            fps = np.stack(fps) # NxN_BITS
            freqs = fps.sum(0)/fps.shape[0]
            fps = freqs > conserved_common_bits_fraction
            assert np.count_nonzero(fps) > 0, "Error, query inputs %s have not conserved bits at fraction > %f "%(
                str(qsmi_list), conserved_common_bits_fraction)
        else:
            fps = np.concatenate(fps) #For Tversky aginst several smiles, we have NxNUM_BITS bits
        return ",".join(qsmi_list), fps

    query_fps = db.from_sequence(query_smi_list).map(compute_query_fp).filter(None.__ne__)

    query_smi_list, query_fps,  = zip(*query_fps.compute())
    query_smi_list = list(query_smi_list)
    query_fps = list(query_fps)

    query_fps = np.stack(query_fps, axis=0) # N_queriesxfingerprint_n_bits type boolean (wasting 7/8 memory but faster)
    resultsDict, n_found = search_fps_matrix(query_fps, database_dir, n_hits_per_smi=n_hits_per_smi,
                                             query_smi_list=query_smi_list, output_name=output_name, backend=backend,
                                             metric=metric, verbose=verbose, n_cpus=n_cpus)
    total_time = time.time() - starting_time
    if verbose: print("Total time for %d smi: %s" % ( n_found, str(datetime.timedelta(seconds=total_time))))
    return resultsDict

def search_fps_matrix(query_fps, database_dir, n_hits_per_smi=30, query_smi_list=None, output_name=None, backend="numpy",
                    metric="Tanimoto", verbose=True, n_cpus=1):

    matched_similarities= np.ones( (len(query_fps), n_hits_per_smi) ) * -np.inf              #query_id, hit_num, similarity
    matched_ids = np.ones( (len(query_fps), n_hits_per_smi, 2), dtype= np.int64 ) * -1  #query_id, hit_num, [ file_id, hit_id]
    if query_smi_list is None:
        query_smi_list = ["q%d"%i for i in range(query_fps.shape[0])]
    kwargs = {}

    compounds_dbname = os.path.join(database_dir, "compounds.sqlite")

    if backend == "numpy":
         process_chunk_fun = process_chunk_using_numpy
         assert metric == "Tanimoto", "Error, numpy backend only supports Tanimoto metric"
         assert isinstance(query_smi_list[0],str) or len(query_smi_list[0]) ==1, "Error, numpy backend is not available when computing similarity with multiplicity > 1"
    elif backend == "numba":
        process_chunk_fun = process_chunk_using_numba

        def get_db_fp_weights():
            con = sqlite3.connect(compounds_dbname)
            bit_counts = pd.read_sql_query("SELECT bit_counts FROM bit_counts", con)["bit_counts"].values
            con.close()
            bit_counts += 1
            bit_logprob = np.log(bit_counts/bit_counts.sum())
            return bit_logprob

        def get_db_fp_logProb2D():
            con = sqlite3.connect(compounds_dbname)
            cur = con.cursor()
            cur.execute("SELECT bit_counts_2d_csv FROM bit_counts_2d_csv")
            csv_bits = cur.fetchone()[0]
            bit_counts2D = pd.read_csv(io.StringIO(csv_bits)).values
            con.close()

            bit_counts2D += 1
            p_a_and_b = bit_counts2D/bit_counts2D.sum()
            p_a = np.diag(bit_counts2D)
            p_a = p_a / p_a.sum()
            p_a_if_b = np.log(p_a_and_b/p_a) #P(a|b)
            return p_a_if_b

        if metric == "Tanimoto":
            kwargs["metric_fun"] = jaccard_numba
        elif metric == "TanimotoW":
            bit_log_weights = 1./get_db_fp_weights()
            @numba.jit(**config.NUMBA_kwARGS)
            def _jaccard_weighted_numba(query_fp, db_fp):
                return jaccard_weighted_numba(query_fp, db_fp, bit_log_weights)
            kwargs["metric_fun"] = _jaccard_weighted_numba

        elif metric == "Tversky":
            kwargs["metric_fun"] = tversky_numba
        elif metric == "Fp_bits_frequency":
            bit_logprob = get_db_fp_weights()
            @numba.jit(**config.NUMBA_kwARGS)
            def _fp_bits_frequency(query_fp, db_fp):
                return fp_bits_frequency(query_fp, db_fp, bit_logprob)
            kwargs["metric_fun"] = _fp_bits_frequency
        elif metric == "Fp_bits_frequency2D":
            bit_log_prob2D = get_db_fp_logProb2D()
            @numba.jit(**config.NUMBA_kwARGS)
            def _fp_bits_frequency2D(query_fp, db_fp):
                return fp_bits_frequency2D(query_fp, db_fp, bit_log_prob2D)
            kwargs["metric_fun"] =_fp_bits_frequency2D

        elif callable(metric):
            kwargs["metric_fun"] = metric
        else:
            raise ValueError(f"Metric {metric} is not a valid option")

    else:
        raise ValueError("Error, not valid backend %s"%backend)

    def process_one_subFile(fname):
        return process_one_subFile_by_chunks(query_fps, fname, n_hits_per_smi,
                                              process_chunk_fun=process_chunk_fun,
                                              verbose=verbose, **kwargs)

    fingerprints_dir = os.path.join(database_dir, "fingerprints")
    filenames = filter( lambda x:  x.endswith(".fingerprints.BitVect"), sorted(os.listdir(fingerprints_dir)))
    filenames = list( map( lambda x:  os.path.join(fingerprints_dir,x), filenames ) )

    if config.USE_DASK_FOR_SEARCH:
        bag = db.from_sequence(enumerate(filenames)).map(process_one_subFile).fold(combine_two_chunk_searches,
                                                                    initial=(matched_similarities, matched_ids))
        matched_similarities, matched_ids = bag.compute()
    else:
        from joblib import Parallel, delayed
        from functools import reduce
        from joblib.externals.loky import set_loky_pickler
        from joblib import wrap_non_picklable_objects
        process_one_subFile = wrap_non_picklable_objects(process_one_subFile)
        set_loky_pickler('dill')
        all_subpartitions = Parallel(n_jobs=n_cpus, backend="loky", verbose=verbose)(delayed(process_one_subFile)( (i, fname) ) for i,fname in enumerate(filenames))
        matched_similarities, matched_ids = reduce(combine_two_chunk_searches, all_subpartitions, (matched_similarities, matched_ids) )
    if verbose: print("Binary search completed! Looking into database")
    # print( matched_similarities )
    # print( matched_ids )



    basenames = [ os.path.basename(fname).split(".")[0] for fname in filenames]
    # matches_db_entries_compound = [ (int(rowNum), basenames[fileNum]) for fileNum, rowNum in matched_ids.reshape(-1, 2) ]

    con = sqlite3.connect(compounds_dbname)
    cur = con.cursor()

    resultsDict = OrderedDict([])
    n_found = 0
    for queryIdx in range(matched_ids.shape[0]):
        matches_db_entries_compound = [(int(rowNum), basenames[fileNum]) for fileNum, rowNum in matched_ids[queryIdx,...].reshape(-1, 2) if rowNum>=0]
        matches_sim_compound = matched_similarities[queryIdx,...]
        query_smi = query_smi_list[queryIdx]
        if verbose: print("> Query: %s"%query_smi)
        # print( matches_db_entries_compound )
        if len(matches_db_entries_compound) > 0:
            resultsDict[query_smi] = []
            found_in_db=0
            repeated = set([])
            for sim, match_entry in sorted(zip(matches_sim_compound, matches_db_entries_compound), key= lambda x: -x[0]):
                # print( match_entry )
                for row in cur.execute("SELECT compoundId, smi FROM smiles NATURAL JOIN ("
                                       " SELECT compoundId FROM compounds WHERE rowNum=? AND fileSource=?)", match_entry):
                    cid_smi = (row[0], row[1])
                    if cid_smi in repeated:
                      continue
                    repeated.add(cid_smi)
                    resultsDict[query_smi].append((sim, row[0], row[1]))
                    if verbose: print("%.4f\t%s\t%s" % (sim, row[0], row[1]))
                    n_found +=1
                    found_in_db+=1
            if len(repeated) == len(matches_sim_compound):
                assert found_in_db == len(matches_db_entries_compound), "Error in database!!"
        else:
            if verbose: print("No matching compounds were found")

    con.close()
    for key in resultsDict:
        resultsDict[key].sort(key=lambda x: -x[0])

    if output_name:
        with open(output_name, "w") as f:
            json.dump(resultsDict, f)  # >>> data = json.loads(json_str, object_pairs_hook=OrderedDict)
    return resultsDict, n_found

def add_common_arguments(parser):
    import sys
    from argparse import FileType

    parser.add_argument('smiles_query', type=FileType('r'), default=sys.stdin,  # nargs=None,
                        help="smiles file with as many smiles as rows. For Tversky similarity, two smiles per row,"
                             " comma separated, can be provided. For combined fingerprint search ('Fp_bits_frequency'),"
                             " separate as many smiles as you want with commas")

    parser.add_argument('-m', '--metric', choices=["Tanimoto", "Tversky", "Fp_bits_frequency"], default="Tanimoto",
                        required=False, help="metric to use")

    parser.add_argument('-c', '--combine_fingerprints', action="store_true", default=False,
                        help="Combine all fingerprints of the smiles into one single fingerprint ")

    parser.add_argument('--conserved_common_bits_fraction', type=float, default=0.5,
                        help="Select as conserved bits for combine_fingerprints those shared by at least %(default)"
                             " fraction of smiles ")

    parser.add_argument('-n', '--n_hits_per_smi', type=int, default=30,
                        help="K highest score molecules to retrieve per query smiles ")

    parser.add_argument('-b', '--backend', choices=["numpy", "numba"], default="numba", required=False,
                        help="computation backend ")

    parser.add_argument('--n_cpus', type=int, default=1,  help="number of cpus to use ")

    parser.add_argument('-o', '--output_name', type=str, required=True,
                        help="The fname for a json file where search results will be stored")

    parser.add_argument('-v', '--verbose', action="store_true", default=False,
                        help="Print to stdout working information ")

def mainSearch():
    from argparse import ArgumentParser
    parser = ArgumentParser(prog="fast_similarity_search", description="Find the K most similar compounds in the database. ")
    add_common_arguments(parser)

    parser.add_argument('-d', '--database_dir', help="the directory where compounds database was compiled",  required=True)
    args = parser.parse_args()
    query_smi_list =  args.smiles_query.read().splitlines()
    query_smi_list = [ x.strip().split(",") for x in query_smi_list]
    assert len(query_smi_list) >0, "Error, no smiles provided"
    if args.combine_fingerprints:
        assert 0 < args.conserved_common_bits_fraction < 1, "Error 0 < conserved_common_bits_fraction < 1 required"
    # print(query_smi_list)
    if config.USE_DASK_FOR_SEARCH:
        dask_client = get_parallel_client()
        print(dask_client.dashboard_link)

    # numba.set_num_threads(multiprocessing.cpu_count()/n_dask_workers)
    search_smi_list(query_smi_list, args.database_dir, n_hits_per_smi=args.n_hits_per_smi, output_name=args.output_name,
                    conserved_common_bits_fraction=args.conserved_common_bits_fraction,
                    backend=args.backend, metric = args.metric, n_cpus=args.n_cpus,
                    verbose=args.verbose)
    if config.USE_DASK_FOR_SEARCH:
        del dask_client
    print("%s search done!"%args.database_dir)



if __name__ == "__main__":
    mainSearch()

    '''

echo -e "CC1CCSCCN1CCC1=CN=CC(F)=C1\nCOC\nC1C(C(C(C(O1)O)O)O)O" | python -m similaritySearch.similarity_searcher_search_onePartition -d ~/oxford/enamine/fingerprints_db --n_cpus 2 -o ~/tmp/simiSearch/first_partition.json -

tail  -n 4 ~/oxford/enamine/cxsmiles/big_example1.txt | cut -f 1 | python -m similaritySearch.similarity_searcher_search_onePartition -d ~/oxford/enamine/fingerprints_db --n_cpus 2 -o ~/tmp/simiSearch/first_partition.json -v -

echo -e "CC1CCSCCN1CCC1=CN=CC(F)=C1,C1CCSCCN1CCC1=CN=CC(Cl)=C1\nCOCCCCCN,C1C(C(C(C(O1)O)O)O)O" | python -m similaritySearch.similarity_searcher_search_onePartition -d ~/oxford/enamine/fingerprints_db --n_cpus 2 -o ~/tmp/simiSearch/first_partition.json -


    '''
