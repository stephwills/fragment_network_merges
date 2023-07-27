#!/usr/bin/env python3
import argparse
import os
import re
import tempfile
import warnings
from subprocess import check_output

import fragment_network_merges.similaritySearch.similaritySearchConfig as config

BASH_TEMPLATE='''###################
%(conda_activate)s
%(env_vars)s
%(cmd)s

###################
'''



CONDOR_TEMPLATE='''###################
Executable      = /usr/bin/bash
Arguments       = %(bash_tmpFile)s
Universe        = vanilla
Output          = %(logdirs)s/$(Cluster).$(Process).out
Error           = %(logdirs)s/$(Cluster).$(Process).err
Log             = %(logdirs)s/$(Cluster).$(Process).log

request_cpus   = %(n_cpus)s
%(additional_requirements)s

Queue

###################

'''


'''
#additional_requirements example:

#request_memory = 4096
#request_disk   = 262144
#request_gpus    = 1

#requirements = (TARGET.Machine == "pulsar-exec-node-gpu-2.xchem.novalocal")

###################
'''

def submit_to_condor(cmd, n_cpus, memory=None, gpus=None, nodename=None, env_vars=None,
         logdirs=config.DEFAULT_LOGS_DIR, tmpdir=config.DEFAULT_SUBMITS_DIR, conda_activate=None,
         only_print=False, **kwargs):
    if len(kwargs) >0:
        warnings.warn("Some kwargs were not used %s"%str(kwargs))

    args = dict(cmd=cmd, n_cpus=n_cpus, tmpdir=tmpdir, logdirs=logdirs)

    if env_vars:
        env_vars = ""
        for env_var in args["env_vars"]:
            env_vars+= ("export " +env_var+ "\n")
    else:
        env_vars = ""
    args["env_vars"]= env_vars

    if conda_activate:
        args["conda_activate"] ='''
eval "$(conda shell.bash hook)"
conda activate %s
        '''%conda_activate
    else:
        args["conda_activate"] = ""

    bash_str = BASH_TEMPLATE%(args)

    additional_requirements=""
    if memory:
        additional_requirements += "request_memory = %d\n"%memory
    if gpus:
        additional_requirements += "request_gpus = %d\n"%gpus
    if nodename:
        additional_requirements += 'requirements = (TARGET.Machine == "%s")\n'%nodename

    args["additional_requirements"] = additional_requirements

    if only_print:
        print("***********************************")
        print(bash_str )
        print("***********************************")

    with tempfile.TemporaryDirectory(dir=tmpdir) as tmp:
        with tempfile.NamedTemporaryFile(mode="w", dir=tmpdir, suffix="_launch.sh", delete=False) as f:
            bash_tmpFile = f.name
            f.write( bash_str)
        args["bash_tmpFile"] = bash_tmpFile
        condor_str = CONDOR_TEMPLATE%(args)

        if only_print:
            print("***********************************")
            print(condor_str)
            print("***********************************")
        else:
            condor_tmpFile = os.path.join(tmp, "launch.condor")
            with open( condor_tmpFile, "w") as tmpfile:
                tmpfile.write( condor_str)

            out = check_output(["condor_submit", condor_tmpFile ])
            match = re.match(re.compile(".*cluster (\d+).*", re.DOTALL), out.decode())
            if match:
                jobId = match.group(1)
            else:
                raise Exception("Error submitting job")
            return jobId

if __name__ == "__main__":

    parser = argparse.ArgumentParser("utility to send commands to condor queue")

    parser.add_argument("--n_cpus", type=int, required=True, help="number of cpus")
    parser.add_argument("--memory", type=int, required=False, default=None, help="Total memory in MB. Default %(default)s")
    parser.add_argument("--gpus", type=int, required=False, default=None, help="Number of GPUs")
    parser.add_argument("--nodename", type=str, required=False, default=None, help="node where job will be executed")

    # parser.add_argument("--bindir", type=str, required=False, default=None, help="directory where the binary lives") #TODO
    parser.add_argument("--logdirs", type=str, required=False, default=config.DEFAULT_LOGS_DIR,
                        help="Logs directory. Default %(default)s")
    parser.add_argument("--tmpdir", type=str, required=False, default=config.DEFAULT_SUBMITS_DIR,
                        help="Logs directory. Default %(default)s")

    parser.add_argument("--env_vars", type=str, nargs="+", required=False, default=[], help="enviramental variables")
    parser.add_argument("--conda_activate", type=str, required=False, default="", help="the name of the environment to activate") #TODO: this is not working properly

    parser.add_argument("--only_print", action="store_true", help="print files instead submitting")

    parser.add_argument("cmd", type=str, help="commands between \"\"")

    args = vars(parser.parse_args())
    print( args )
    jobId = submit_to_condor(**args)
    print(jobId)

'''
python -m utils.send_to_condor --n_cpus 4 "python -c 'print(0)'"
'''