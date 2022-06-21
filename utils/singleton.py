"""Singleton metaclass"""


class Singleton(type):
    """
    Metaclass to prevent several instances of the class Config
    """
    _instances = {}

    def __call__(cls, *args, **kwargs):
        if cls not in cls._instances:
            cls._instances[cls] = super(Singleton, cls).__call__(*args, **kwargs)
        return cls._instances[cls]
