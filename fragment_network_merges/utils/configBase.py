import os
from fragment_network_merges.utils.singleton import Singleton


class ConfigBase(metaclass=Singleton):
    CONFIG_DICT = dict()

    def get(self, key):
        val = os.environ.get(key, None)
        if val:
            val = self.default_types[key](val)
            return val
        else:
            return self.CONFIG_DICT[key]

    # this method will allow instances to access config properties as obj.PROPERTY; useful for subclassing
    def __getattr__(self, key):
        if key in self.__dict__:
            val = self.__dict__[key]
            return val
        else:
            return self.get(key)

    def __setattr__(self, key, value):
        if value is not None:
            super(ConfigBase, self).__setattr__(key, value)
        # self.CONFIG_DICT[key] = value

    def __init__(self):
        self.default_types = {}
        for k, v in self.CONFIG_DICT.items():
            if v is None:
                _type=str
            else:
                _type = type(v)
            self.default_types[k] = _type


config_base = ConfigBase()
