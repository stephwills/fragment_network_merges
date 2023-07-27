import os

from fragment_network_merges.utils.configBase import ConfigBase

class Config(ConfigBase):
    CONFIG_DICT = dict(
        FRAGALYSIS_DATA_DIR="/home/swills/Oxford/data/Fragalysis",
        WORKING_DIR=os.path.join(os.getcwd(), "data"),
        N_CPUS_FILTER_PAIR=1,
        USE_NEO4J_INSTEAD_API=True,
        # ONLY FOR USE_NEO4J_INSTEAD_API=True
        NEO4J_URI="bolt://localhost:7687",
        NEO4J_USER="swills",
        NEO4J_PASS=None,  # export to bash before execution; export NEO4J_PASS="myPass",
        # ONLY FOR USE_NEO4J_INSTEAD_API=False
        RESTAPI_USER="swills",
        RESTAPI_PASS=None,  # export to bash before execution; export RESTAPI_PASS="myPass",
        # FOR NEO4J DATABASE QUERY
        NUM_HOPS=2,  # number of hops specified in the merging query,
        MIN_HAC=15,  # min number of heavy atoms of merges to retrieve
        MAX_HAC=1000,  # max number of heavy atoms of merges to retrieve
        MIN_HAC_FA=7,  # min number of heavy atoms of fragment A before expansion
        # FOR FILTERING THE FRAGMENTS BEFORE QUERYING THE DATABASE
        MAX_FRAG_DIST=5,  # the maximum distance between two fragments for them to be merges (angstroms)
        # FOR FILTERING THE SYNTHONS BEFORE QUERYING THE DATABASE
        MIN_CARBONS=3,  # the minimum number of carbons for a synthon to be used in the expansion
        SYNTH_OVERLAP=0.5,  # the proportion of overlap between a synthon found in both fragments for it to be ruled out
        NUM_ATOMS_NOT_MCS=3,  # number of atoms not included in MCS between fragments when checking for too similar pairs
        MCS_RMSD_THRESHOLD=2,  # threshold for ruling out fragment pairs that are too similar (RMSD between MCS atoms)
    )

config_merge = Config()
