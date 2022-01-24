import os

from merge.utils import Singleton


class Config(metaclass=Singleton):
    CONFIG_DICT = dict(
        FRAGALYSIS_DATA_DIR="/home/swills/Oxford/data/Fragalysis",
        WORKING_DIR=os.path.join(os.getcwd(), "data"),
        N_CPUS_FILTER_PAIR=2,
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
        MAX_HAC=30,  # max number of heavy atoms of merges to retrieve
        # FOR FILTERING THE FRAGMENTS BEFORE QUERYING THE DATABASE
        MAX_FRAG_DIST = 5, # the maximum distance between two fragments for them to be merges (angstroms)
        # FOR FILTERING THE SYNTHONS BEFORE QUERYING THE DATABASE
        MIN_CARBONS=3,  # the minimum number of carbons for a synthon to be used in the expansion
        SYNTH_OVERLAP=0.5,  # the proportion of overlap between a synthon found in both fragments for it to be ruled out
    )

    @classmethod
    def get(cls, key):
        val = os.environ.get(key, None)
        if val:
            return val
        else:
            return cls.CONFIG_DICT[key]

    # this method will allow instances to access config properties as obj.PROPERTY; useful for subclassing
    def __getattr__(self, key):
        if key in self.__dict__:
            val = self.__dict__[key]
            return val
        else:
            return self.get(key)

    def __setattr__(self, key, value):
        if value is not None:
            super(Config, self).__setattr__(key, value)
        # self.CONFIG_DICT[key] = value
        
    def __init__(self):
        pass


config_merge = Config()
