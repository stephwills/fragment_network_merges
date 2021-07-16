# Fragment Network merges

## Motivation

Fragments can be optimized to become lead-like compounds by fragment merging, whereby fragments that bind to a target in adjoining or partially overlapping space are merged by finding compounds that incorporate structural motifs from each. The Fragment Network is a graph database containing catalogue compounds. Given there are two fragment hits that exist as nodes in the network, the database can be queried to find catalogue compounds that represent merges of the two.

## Features

fragment_network_merges provides a pipeline for querying the Fragment Network to find fragment merges and filtering the results to find the most promising compounds. The merges require filtering to find those in which the binding mode and interactions made by the original fragments are most likely to be conserved. The filtering steps implemented here are as follows:

* **Descriptor filter:** Merges are filtered using Lipinski’s rules.
* **Expansion filter:** Merges that represent expansions (i.e. elaborations of a single fragment rather than a true merge) are removed
* **Embedding filter:** Constrained embedding is performed whereby a conformation is generated for the merge in which the atoms from the parent fragments are constrained using their original coordinates; merges are removed if it is impossible to generate a physically reasonable structure
* **Overlap filter:** Merges are removed if they do not fit the protein pocket (i.e. >15% of the volume of the merge ‘overlaps’ with the protein structure
* **Fragmenstein filter:** Merges are ‘placed’ in the protein using Fragmenstein, which uses the MCS between the merge and the original fragments and performs minimization using PyRosetta
  * Merges are removed if both fragments are not able to be used for placement, the combined RMSD is >1Å and ΔΔG is positive
* **Scoring:** Merges are scored according to the proportion of interactions made by the original fragments that are maintained by the merge (calculated using ODDT’s PLIF implementation); other scoring metrics can also be used for further filtering and prioritization (e.g. SuCOS)
