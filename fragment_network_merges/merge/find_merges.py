"""
Returns the appropriate fragment network searcher depending on whether port forwarding or RestAPI used.
"""

from fragment_network_merges.merge.config_merge import config_merge
from fragment_network_merges.merge.find_merges_neo4j import MergerFinder_neo4j
from fragment_network_merges.merge.find_merges_restAPI import MergerFinder_restAPI


def getFragmentNetworkSearcher(
    use_neo4j_instead_api=config_merge.USE_NEO4J_INSTEAD_API, **kwargs
):
    """
    Returns the Fragment Network searcher according to whether neo4j or restAPI are used.
    """
    if use_neo4j_instead_api:
        return MergerFinder_neo4j(**kwargs)
    else:
        return MergerFinder_restAPI(**kwargs)
