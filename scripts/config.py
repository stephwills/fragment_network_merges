import os

class Config():
    CONFIG_DICT = dict(
        FRAGALYSIS_DATA_DIR = os.path.expanduser("~/oxford/data/fragalysis"),
        WORKING_DIR = os.path.join(os.getcwd(), "data") # I would prefer this as part of the argparse
    )
    
    @classmethod
    def get(cls, key):
        val = os.environ.get(key, None)
        if val:
            return val
        else:
            return cls.CONFIG_DICT[key]

