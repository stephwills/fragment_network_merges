# Fragment Network merges

## Motivation

Fragments can be optimized to become lead-like compounds by fragment merging, whereby fragments that bind to a target in adjoining or partially overlapping space are merged by finding compounds that incorporate structural motifs from each. The Fragment Network is a graph database containing catalogue compounds. Given there are two fragment hits that exist as nodes in the network, the database can be queried to find purchasable compounds that represent merges.

## Overview

fragment_network_merges provides a pipeline for querying the Fragment Network to find fragment merges and filtering the results to find the most promising compounds.
There are two packages, `merge` and `filter`. `merge` contains scripts to query the neo4j database and retrieve merges. `filter` contains scripts that filter the merges down to a manageable number of the most promising compounds.

### Querying
Given there are two fragments to merge, the database is queried using one of the fragments as the seed fragment. All nodes within a specified number of hops (representing close neighbours) of the seed fragment are identified (the default number of hops is 2). The query then identifies expansions of these compounds that incorporate a substructure from the other fragment in the pair. All possible substructures are enumerated and are queried for in the database. The querying process is asymmetric; thus, there are two queries for each pair of fragments. 

### Filtering pipeline
The merges require filtering to find the most promising compounds and those in which the binding mode and interactions made by the original fragments are most likely to be conserved. The filters have been implemented using an abstract class and thus new filters can be added to the pipeline. The pipeline also has the flexibility to allow the user to choose which filters to run and in what order. There are also scoring functions implemented that can be used to score filtered compounds and allow further prioritization in future steps. The available filtering and scoring functions are as follows:

* **Descriptor filter:** Merges are filtered using Lipinski’s rules and using a rotatable bond count limit of 10.
* **Non-ring bond filter:** Filters out 'stringy' molecules with long linkers depending on the number of consecutive non-ring bonds.
* **Expansion filter:** Merges that represent expansions (i.e. elaborations of a single fragment rather than a true merge) are removed.
* **Embedding filter:** Constrained embedding is performed whereby a conformation is generated for the merge in which the atoms from the parent fragments are constrained using their original coordinates; merges are removed if it is impossible to generate a physically reasonable structure.
* **Energy filter:** Filters out compounds with unrealistic conformations depending on the ratio between the energy of the constrained conformation with the energy of several unconstrained conformations.
* **Overlap filter:** Merges are removed if they do not fit the protein pocket, determined by calculating the clash distance between the ligand and the protein. 
* **Fragmenstein filter:** Merges are ‘placed’ in the protein using [Fragmenstein](https://github.com/matteoferla/Fragmenstein), which uses the MCS between the merge and the original fragments and performs minimization using PyRosetta.
  * Merges are removed if both fragments are not able to be used for placement, the combined RMSD is >1Å and ΔΔG is positive.
* **PLIP interaction score:** Calculates interactions using the Protein-Ligand Interaction Profiler (PLIP) and scores them according to the proportion of interactions made by the original fragments that are maintained by the merge.
* **SuCOS score:** Calculates the mean shape and pharmacophoric overlap with the original fragments.

## Installation

Create environment from the environment.yml file using `conda env create -f environment.yml`

Clone the GitHub repository using `git clone https://github.com/stephwills/fragment_network_merges.git` and install using `pip install -e .`

### Optional

Packages for optional filtering steps include:
* plip == 2.2.2
* fragmenstein == 0.9.11

To install PLIP, follow the instructions on the GitHub repository (requires OpenBabel). To install Fragmenstein, follow the instructions in the GitHub repository. Fragmenstein requires PyRosetta, which requires an academic license.
## Usage

### Merge

Querying the database can be done by port-forwarding with Kubernetes. The configuration file `merge/config_merge.py` should be edited before running for the first time.

To use neo4j with Kubernetes, you need to access the network-db using port forwarding. In the config file, set USE_NEO4J_INSTEAD_API to True or False depending on which option is being used. The neo4j username should also be set (NEO4J_USER or RESTAPI_USER). The user's neo4j password should be exported to bash before running (i.e. `export NEO4J_PASS="myPass"`).

At this moment port forwarding can be executed as follows:

```
export KUBECONFIG=$HOME/.kube/_fragment_network_config
kubectl port-forward pods/graph-0 7474:7474  &
kubectl port-forward pods/graph-0 7687:7687 &
```

In order to define the directory pointing to the data, set FRAGALYSIS_DATA_DIR in `merge/config_merge.py`. The FRAGALYSIS_DATA_DIR should contain the target data in the Fragalysis format.

Several parameters for querying can also be set in the config file. These include parameters involved in the querying the database (e.g. number of optional hops to be made), parameters involved in filtering the fragment pairs and parameters involved in filtering the synthons used for expansion.

To run the query, `merge/query.py` can be run from the command line, specifying the target, the list of fragments to be merged and the directory to save the files. The options are shown below:

```
usage: query.py [-h] -t TARGET -f FRAGMENTS [FRAGMENTS ...] -o OUTPUT_DIR
                [-w WORKING_DIR]

optional arguments:
  -h, --help            show this help message and exit
  -t TARGET, --target TARGET
                        the protein target (e.g. nsp13)
  -f FRAGMENTS [FRAGMENTS ...], --fragments FRAGMENTS [FRAGMENTS ...]
                        the list of fragments to merge. E.g. x0032_0A x0034_0B
  -o OUTPUT_DIR, --output_dir OUTPUT_DIR
                        the directory to write the merge files to
  -w WORKING_DIR, --working_dir WORKING_DIR
                        the directory where intermediate results will be
                        computed

Example python -m merge.query -t nsp13 -f x0034_0B x0176_0B x0212_0B -o
data/example_folder
```
This will enumerate all possible pairs of fragments and for each pair, will query the database to find the substructures of one fragment for expansion, and then find expansions of the seed fragment using each substructure. For each pair, a json file is created in the specified output directory, named according to the merge, e.g. `x0034_0B_x0176_0B.json` will include all merges found involving expansions of fragment 34 using substructures of fragment 176. These are saved in a dictionary, with the dictionary keys indicating the substructure (also referred to as 'synthon') used in the expansion, and the values representing the list of merges found (in SMILES format).

### Filter

**Filtering merges**

The config file `filter/config_filter.py` should be edited prior to running the filtering pipeline.
In order to define the directory pointing to the data, set FRAGALYSIS_DATA_DIR in `filter/config_filter.py`. The FRAGALYSIS_DATA_DIR should contain the target data in the Fragalysis format.

Merges can be passed through a series of filters to reduce them to a more manageable number. Individual filters exist in separate modules and are built upon the abstract Generic_Filter class (in `generic_filter.py`). 
The pipeline is flexible to allow the filters to be used in a chosen order, for filters to be removed or for new filters to be implemented. However, the filters should follow a logical order; i.e. if a filter requires a 3D conformation of the molecule then there must be a filter prior in which a conformation has been generated.
The filtering pipeline is specified by FILTER_PIPELINE in the config file; the class names are entered into a list in the order in which the user wants to run the filters.

**Scoring merges**

Filtered molecules can then be scored using scoring functions built upon the abstract class Score_generic in `generic_scoring.py`. Thus far, an interaction score (`plip_ifp_score.py`) and SuCOS score (`sucos_score.py`) are available. 
The scoring functions to be used are specified by SCORING_PIPELINE in the config file.

Other thresholds used in the filters can also be specified in the config file.

The filtering pipeline can be run from the command line using `filter_pipeline.py`. The pipeline is run for a single fragment 'pair' (remembering that the expansion technique is asymmetrical). The pipeline takes as input an individual json file produced by the merging pipeline above.

The options are shown below:

```
usage: filter_pipeline.py [-h] [-f MERGE_FILE] [-a FRAGMENTA] [-b FRAGMENTB]
                          [-t TARGET] [-o OUTPUT_DIR] [-w WORKING_DIR] [-s]

optional arguments:
  -h, --help            show this help message and exit
  -f MERGE_FILE, --merge_file MERGE_FILE
                        the json file containing the merges
  -a FRAGMENTA, --fragmentA FRAGMENTA
                        the name of fragment A, e.g. x0107_0A
  -b FRAGMENTB, --fragmentB FRAGMENTB
                        the name of fragment B, e.g. x0434_0A
  -t TARGET, --target TARGET
                        the name of the target, e.g. nsp13 or Mpro
  -o OUTPUT_DIR, --output_dir OUTPUT_DIR
                        the directory to write the filtered files to
  -w WORKING_DIR, --working_dir WORKING_DIR
                        the directory to write the intermediate files to
  -s, --sim_search      whether the data is similarity search data or not

Example python -m filter.filter_pipeline -f x0034_0B_x0176_0B.json -m
x0034_0B_x0176_0B -t nsp13 -o data/example_folder
```

Within the working directory a temporary sub-directory is created, where all the Fragmenstein files are temporarily written to. The files of merges that pass the filter are moved to the output directory.

Running the filtering script creates two json files, a `MERGE_filtered.json` and a `MERGE_failures.json`. The former contains a dictionary containing the successfully filtered compounds, their unique name (which is used to name the Fragmenstein files), SMILES, substructure used for expansion and any scores used. The latter file contains the compounds that did not pass the filters, their unique name, the SMILES and which filter they failed at. Other files containing information on the timings for running each filter and Fragmenstein timings are also created.

## Tests

Unit tests are available for testing the individual filtering and scoring functions. 