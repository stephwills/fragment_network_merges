import os

from scripts.utils import Singleton


class Config(metaclass=Singleton):
    CONFIG_DICT = dict(
        FRAGALYSIS_DATA_DIR = os.path.expanduser("~/oxford/data/fragalysis"),
        WORKING_DIR = os.path.join(os.getcwd(), "data"), # I would prefer this as part of the argparse, but it is better than nothing
        N_CPUS_FILTER_PAIR=2,
        USE_NEO4J_INSTEAD_API=True,
        NEO4J_URI="bolt://localhost:7687",
        NEO4J_USER ="rgarcia",
        NEO4J_PASS=None, #This should never be plain text. Better exporting to bash before execution. export NEO4J_PASS="myPass",
        RESTAPI_USER = "rgarcia",
        RESTAPI_PASS = None, #This should never be plain text. Better exporting to bash before execution. export NEO4J_PASS="myPass"
    )
    
    @classmethod
    def get(cls, key):
        val = os.environ.get(key, None)
        if val:
            return val
        else:
            return cls.CONFIG_DICT[key]

    #This method will allow instances to access config properties as obj.PROPERTY. Useful for subclassing
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
        


config = Config()

if __name__ == "__main__":

    conf = Config()

    print(conf.N_CPUS_FILTER_PAIR)
    conf.N_CPUS_FILTER_PAIR = 3
    print(conf.N_CPUS_FILTER_PAIR)
    del conf.N_CPUS_FILTER_PAIR
    print(conf.N_CPUS_FILTER_PAIR)
    
    
