"""
Used to generate merges between two fragments using the fragment network.
Uses code from https://github.com/tdudgeon/fragment-network-merges.
"""

from scripts.config import config
from scripts.find_merges_neo4j import MergerFinder_neo4j
from scripts.find_merges_restAPI import MergerFinder_restAPI


def getFragmentNetworkSearcher(use_neo4j_instead_api= config.USE_NEO4J_INSTEAD_API, **kwargs):
    if use_neo4j_instead_api:
        return MergerFinder_neo4j(**kwargs)
    else:
        return MergerFinder_restAPI(**kwargs)
