#TODO: This config file is a quick and dirty solution. It requires to be integrated using standard config /argument parsers
import os
import tempfile

FINGERPRINT_TYPE="morgan" #pharmacophore
FINGERPRINT_NBITS=2048
FINGERPRINT_RADIUS=3
FINGERPRINT_USE_CHIRALITY=False
FINGERPRINT_USE_FEATURES=True
FINGERPRINT_USE_BOND_TYPES=True


N_LINES_PER_CHUNK = int(5e5) # int(1e6)
N_CPUS = 1
DASK_WORKER_MEMORY = '8GB'
TMP_DIR = tempfile.gettempdir()
DISABLE_DASHBOARD=True
USE_DASK_FOR_SEARCH= False

CORRECT_METRIC_BY_BITS_ON_IN_CANDIDATE_FREQUENCY_BASED_METRICS = NotImplemented #False

DASK_DISTRIBUTED__COMM__TIMEOUTS__CONNECT = 100

NUMBA_kwARGS = dict(nopython=True, cache=True, nogil=True, parallel=True) #Warning, if using parallel=True, set os.environ['OPENBLAS_NUM_THREADS'] = '1'

########## Cluster params
WAIT_TIME_LAUNCH_QUEUE = 5

#TODO: Move this to another config file
SHARED_TMP = "/data/xchem-fragalysis/sanchezg/tmp/"
DEFAULT_LOGS_DIR = "/data/xchem-fragalysis/sanchezg/logs/"
DEFAULT_SUBMITS_DIR = "/data/xchem-fragalysis/sanchezg/submits/"
PATH="/data/xchem-fragalysis/sanchezg/app/miniconda3_2/envs/Fragmenstein/bin:$PATH"
PYTHON_CLUSTER="/data/xchem-fragalysis/sanchezg/app/miniconda3_2/envs/Fragmenstein/bin/python"
