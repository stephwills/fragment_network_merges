#TODO: This config file is a quick and dirty solution. It requires to be integrated using standard config /argument parsers
import os
import tempfile

FINGERPRINT_NBITS=2048
N_LINES_PER_CHUNK = int(5e5) # int(1e6)
N_CPUS = 1
DASK_WORKER_MEMORY = '8GB'
TMP_DIR = tempfile.gettempdir()
DISABLE_DASHBOARD=True
USE_DASK_FOR_SEARCH= False


DASK_DISTRIBUTED__COMM__TIMEOUTS__CONNECT = 100

########## Cluster params
WAIT_TIME_LAUNCH_QUEUE = 5

#TODO: Move this to another config files
SHARED_TMP = "/data/xchem-fragalysis/sanchezg/tmp/"
DEFAULT_LOGS_DIR = "/data/xchem-fragalysis/sanchezg/logs/"
DEFAULT_SUBMITS_DIR = "/data/xchem-fragalysis/sanchezg/submits/"
PATH="/data/xchem-fragalysis/sanchezg/app/miniconda3_2/envs/Fragmenstein/bin:$PATH"
PYTHON_CLUSTER="/data/xchem-fragalysis/sanchezg/app/miniconda3_2/envs/Fragmenstein/bin/python"
