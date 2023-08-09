This describes how to run the `fragment-network-merge-swills23` job from the `comp chem`
category in the `fragment-network` collection.

## What the job does

fragment-network-merge-swills23 provides a pipeline for querying the Fragment Network to find fragment merges and 
filtering the results to find the most promising compound. First, 2D candidates are pulled from a neo4j database.
Then, we filter the merges down to a manageable number of the most promising compounds.

See the [fragment_network_merges documentation](https://github.com/stephwills/fragment_network_merges) for more information.

This work is published in [JCIM](https://doi.org/10.1021/acs.jcim.3c00276).