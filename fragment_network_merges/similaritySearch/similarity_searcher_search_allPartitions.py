import json
import os
import tempfile
import time
from argparse import ArgumentParser
from subprocess import check_call, check_output
import re
import fragment_network_merges.similaritySearch.similaritySearchConfig as config
from fragment_network_merges.similaritySearch.similarity_searcher_collect_results import combine_search_jsons
from fragment_network_merges.utils.send_to_condor import submit_to_condor


def findDB_partitions(db_paths):

  partitions = []
  assert isinstance(db_paths, (list, tuple))
  for db_path in db_paths:
    if os.path.isfile(os.path.join(db_path, "compounds.sqlite")):
      partitions.append(db_path)
    for fname in os.listdir(db_path):
        full_path = os.path.join(db_path, fname)
        if os.path.isfile(os.path.join(full_path, "compounds.sqlite")):
          partitions.append( full_path )
  print("Partitions to search from", partitions)
  return partitions

def parse_memsize(size):
  #https://stackoverflow.com/questions/42865724/parse-human-readable-filesizes-into-bytes/42865957#42865957
  units = {"KB": 2 ** -10, "MB": 1, "GB": 2 ** 10, "TB": 2 ** 20} #UNIT_TO_MBs
  size = size.upper()
  #print("parsing size ", size)
  if not re.match(r' ', size):
      size = re.sub(r'([KMGT]?B)', r' \1', size)
  number, unit = [string.strip() for string in size.split()]
  return int(float(number)*units[unit])

def launch_searcher(run_locally=False, **kwargs):

  partition_name = kwargs["database_dir"]
  output_name = os.path.join( kwargs["working_dir"], os.path.basename(partition_name).split(".")[0]+".json")
  kwargs["output_name"] = output_name
  kwargs["DASK_DISTRIBUTED__COMM__TIMEOUTS__CONNECT"] = config.DASK_DISTRIBUTED__COMM__TIMEOUTS__CONNECT
  kwargs["PATH"] = config.PATH

  if kwargs.get("dask_worker_memory", "-1") == "-1":
    kwargs["dask_worker_memory"] = " "
    kwargs["queue_memory"] = ""
  else:
    kwargs["queue_memory"] = "--memory "+str( parse_memsize(kwargs["dask_worker_memory"])*kwargs["n_cpus"])
    kwargs["dask_worker_memory"] = " DASK_WORKER_MEMORY=%s "%kwargs["dask_worker_memory"]

  cmd_condor = 'python -m utils.send_to_condor --env_vars ' \
        'PATH=%(PATH)s  %(dask_worker_memory)s  ' \
        ' DASK_DISTRIBUTED__COMM__TIMEOUTS__CONNECT=%(DASK_DISTRIBUTED__COMM__TIMEOUTS__CONNECT)s' \
        ' --ncpus %(n_cpus)s %(queue_memory)s  "'

  cmd_args = ' -m similaritySearch.similarity_searcher_search_onePartition ' \
             '-d %(database_dir)s -o %(output_name)s  %(smilesFname)s --n_cpus %(n_cpus)s'

  if "metric" in kwargs:
    cmd_args += " --metric %s "%kwargs["metric"]

  if "n_hits_per_smi" in kwargs:
    cmd_args += " --n_hits_per_smi %s"% kwargs["n_hits_per_smi"]

  if "backend" in kwargs:
    cmd_args += " --backend %s"% kwargs["backend"]

  if kwargs.get("verbose", False):
    cmd_args += " -v "

  if not run_locally:
    python = " " + config.PYTHON_CLUSTER + " "
    # cmd = cmd_condor + python + cmd_args +'"'
    cmd =  python + cmd_args
    cmd = cmd % kwargs
    jobId = submit_to_condor(  cmd , **kwargs )
    print("jobId", jobId)

  else:
    cmd = "%(dask_worker_memory)s " \
          "DASK_DISTRIBUTED__COMM__TIMEOUTS__CONNECT=%(DASK_DISTRIBUTED__COMM__TIMEOUTS__CONNECT)s " \
          "python " + cmd_args

    cmd = cmd%kwargs
    check_call( cmd, shell=True)
    jobId = -1

  print("Running:\n", cmd, "\nis in cluster:", not run_locally)

  return (output_name, jobId)

def globalSearch():
  from fragment_network_merges.similaritySearch.similarity_searcher_search_onePartition import add_common_arguments

  parser = ArgumentParser(prog="fast_similarity_search",
                          description="Find the K most similar compounds in the database")

  add_common_arguments(parser)

  parser.add_argument('-d', '--database_dirs', nargs="+", help="the directory(s) where compounds database was compiled", required=True)

  parser.add_argument('-w', '--working_dir', type=str, required=True,
                      help="The directory where per partition results will be saved")

  parser.add_argument('-l', '--run_locally', action="store_true",
                      help="run computations locally instead submitting to condor") #TODO: change that from boolean to choices {local,condor, slurm...}

  args = parser.parse_args()
  query_smi_str = args.smiles_query.read()
  query_smi_list =query_smi_str.splitlines()
  assert len(query_smi_list) > 0, "Error, no smiles provided"

  kwargs = vars( args )

  if  not kwargs["run_locally"]:
    tmpdir = config.SHARED_TMP
    jobIdAlive = lambda x: True
  else:
    tmpdir = tempfile.tempdir
    def jobIdAlive(jobId):
      out = check_output(["condor_q", "-analyze", str(jobId)])
      print("Is job alive?", out)
      match = re.match(re.compile(".*Job is running.*", re.DOTALL), out.decode())
      if match:
        print("Job %s is alive"%str(jobId))
        return True
      else:
        print("Job %s is not alive"%str(jobId))
        return False

  db_partitions = findDB_partitions(args.database_dirs)
  assert len(db_partitions) > 0, "Error, database partitions not found"

  with tempfile.NamedTemporaryFile(mode="w", dir=tmpdir, suffix=".smi", delete=False) as f:
    f.write( query_smi_str )
    f.flush()
    kwargs["smilesFname"] = f.name

    del kwargs["smiles_query"]
    results_names = []
    for partitionDir in db_partitions:
      kwargs["database_dir"] = partitionDir
      output_name = launch_searcher( **kwargs )
      results_names.append(output_name)

    print("waiting results...")
    while True:
      completed = True
      remaining = 0
      for name, jobId in results_names:
        isCompleted = os.path.isfile(name)
        if not isCompleted and not jobIdAlive(jobId):
          raise Exception("Error, jobId %s is not alive, but results were not computed"%jobId)
        print(name, "jobId", jobId, "is completed: ", isCompleted)
        completed &= isCompleted
        remaining += int(not isCompleted)
      if completed:
        break
      else:
        print("waiting results. Remaining: %d"%(remaining))
        time.sleep(config.WAIT_TIME_LAUNCH_QUEUE)
    output_names, jobIds = zip(*results_names)
    combined_json = combine_search_jsons(output_names)
    with open(args.output_name, "w") as f:
      json.dump(combined_json, f)

if __name__ == "__main__":
  globalSearch()


'''
echo "CCNC(=O)C(F)(F)C(=O)NCC" | python -m similaritySearch.similarity_searcher_search_allPartitions -d ~/oxford/enamine/fingerprints_db ~/oxford/enamine/fingerprints_db2  --run_locally --n_cpus=4 -w ~/tmp/simiSearch/  -o ~/tmp/simiSearch/results.json -

#In pulsar but local
echo "CCNC(=O)C(F)(F)C(=O)NCC" | python -m similaritySearch.similarity_searcher_search_allPartitions -d /data/xchem-fragalysis/sanchezg/oxford/enamine/fingerprints_db /data/xchem-fragalysis/sanchezg/oxford/enamine/fingerprints_db2  --n_cpus=4 -w /data/xchem-fragalysis/sanchezg/tmp/simiSearch/  -o /data/xchem-fragalysis/sanchezg/tmp/simiSearch/results.json -

'''