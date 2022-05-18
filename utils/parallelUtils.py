import os
import re

import dask

from similaritySearch.similaritySearchConfig import N_CPUS, DASK_WORKER_MEMORY, TMP_DIR, DISABLE_DASHBOARD

import os
import logging

from dask.distributed import Client, LocalCluster


# DASK_DISTRIBUTED__COMM__TIMEOUTS__CONNECT="20" #Default is 10

journal = logging.getLogger('Dask_Parallel')
journal.setLevel(logging.DEBUG)

DASK_CLIENT = None

def get_parallel_client(threads_per_worker=None, n_workers=None, memory_limit=None):
    global DASK_CLIENT
    if DASK_CLIENT is None:
        if n_workers is None:
            n_workers = N_CPUS
            threads_per_worker = 1
        if memory_limit is None:
            memory_limit= DASK_WORKER_MEMORY
            if memory_limit == "-1":
                from psutil import virtual_memory
                mem = virtual_memory()
                if mem.total is None:
                    raise ValueError("Error, memory was not determined")
                memory_limit="%dGB"%( (0.9 *mem.total/n_workers) // 2 ** 30)
        dask.config.set({'temporary_directory': os.path.join(TMP_DIR, "dask")})
        if n_workers>1 or threads_per_worker>1:
            if DISABLE_DASHBOARD:
                kwargs = {"dashboard_address":None}
            else:
                kwargs = {}
            DASK_CLIENT = Client( LocalCluster(threads_per_worker=threads_per_worker, n_workers=n_workers, memory_limit=memory_limit, **kwargs)) # dashboard_address=8787
        else:
            DASK_CLIENT = Client( LocalCluster(threads_per_worker=1, n_workers=1) )
        print(DASK_CLIENT, flush=True)
        if not DISABLE_DASHBOARD:
            print(DASK_CLIENT.scheduler_info()['services'])
    return DASK_CLIENT


def apply_func_to_files(folder, file_pattern, function, use_parallel_dask=None, extensions_to_check=None,
                        ids_to_check=None):
    results = []

    if use_parallel_dask is None:
        use_parallel_dask = N_CPUS > 0

    if use_parallel_dask:
        function = dask.delayed(function)
    for root, dirs, fnames in os.walk(folder, followlinks=True):
        for fname in fnames:
            match_obj = re.match(file_pattern, fname)
            if match_obj:
                if extensions_to_check and not os.path.splitext(fname)[-1] in extensions_to_check:
                    raise ValueError(
                        "Error, one of the files found (%s) did not match the required extension (%s)" % (
                            fname, extensions_to_check))
                if ids_to_check:
                    match_id = match_obj.group(1)
                    if match_id not in ids_to_check:
                        continue
                fname = os.path.join(root, fname)
                result = function(fname)
                results.append(result)
    if use_parallel_dask:
        return dask.compute(results)[0]
    else:
        return results
