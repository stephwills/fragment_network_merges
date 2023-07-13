This describes how to run the `fragment_network-merge` jobs from the `comp chem`
category in the `fragment_network` collection.

## What the job does

fragment_network_merges provides a pipeline for querying the Fragment Network to find fragment merges and 
filtering the results to find the most promising compound. First, 2D candidante are pulled from a neo4j database.
Then, we filter the merges down to a manageable number of the most promising compounds.

See the [fragment_network_merges documentation](https://github.com/stephwills/fragment_network_merges) for more information