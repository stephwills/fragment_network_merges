# Fragment Network merges

## Motivation

Fragments can be optimized to become lead-like compounds by fragment merging, whereby fragments that bind to a target in adjoining or partially overlapping space are merged by finding compounds that incorporate structural motifs from each. The Fragment Network is a graph database containing catalogue compounds. Given there are two fragment hits that exist as nodes in the network, the database can be queried to find purchasable compounds that represent merges.

## Pipeline

fragment_network_merges provides a pipeline for querying the Fragment Network to find fragment merges and filtering the results to find the most promising compounds.

### Finding merges
The technique for finding merges is called ‘synthon expansion’, with ‘synthon’ being used here to denote a constituent part of a fragment. If we have two fragments, fragments A and B, this method attempts to find all expansions of fragment A that contain part of fragment B. A series of optional hops can be made before the expansion, in which fragment A gains or loses synthons, to increase diversity in the compounds found. The Fragment Network can be loaded using Neo4j, and the merges are found by querying the database using the cypher query language. 

### Filtering merges
The merges require filtering to find those in which the binding mode and interactions made by the original fragments are most likely to be conserved. The filtering steps implemented here are as follows:

* **Descriptor filter:** Merges are filtered using Lipinski’s rules.
* **Expansion filter:** Merges that represent expansions (i.e. elaborations of a single fragment rather than a true merge) are removed.
* **Embedding filter:** Constrained embedding is performed whereby a conformation is generated for the merge in which the atoms from the parent fragments are constrained using their original coordinates; merges are removed if it is impossible to generate a physically reasonable structure
* **Overlap filter:** Merges are removed if they do not fit the protein pocket (i.e. >15% of the volume of the merge ‘overlaps’ with the protein structure.
* **Fragmenstein filter:** Merges are ‘placed’ in the protein using Fragmenstein, which uses the MCS between the merge and the original fragments and performs minimization using PyRosetta.
  * Merges are removed if both fragments are not able to be used for placement, the combined RMSD is >1Å and ΔΔG is positive.
* **Scoring:** Merges are scored according to the proportion of interactions made by the original fragments that are maintained by the merge (calculated using ODDT’s PLIF implementation); other scoring metrics can also be used for further filtering and prioritization (e.g. SuCOS).

## Installation

To be updated

### Prerequisites

* rdkit
* numpy
* oddt
* neo4j
* biopython 
* Fragmenstein (and PyRosetta)

## Usage

### Querying the database

Querying the database requires access to the Fragment Network with Kubernetes. Prior to running `run_database_query,.py` you need to access the network using port forwarding and to change the username in `find_merges.py` to your own username.

The `preprocessing.py` script requires the data folder containing all the crystal structures (e.g. nsp13, Mpro) to be saved in the fragment_network_merges root directory (this was done to make it easier to access the mol files of the parent fragments). 

To run the query, `run_database_query.py` can be run from the command line, specifying the target, the list of fragments to be merged and the directory to save the files. The options are shown below:

```
usage: run_database_query.py [-h] [-t TARGET] [-f FRAGMENTS [FRAGMENTS ...]]
                             [-o OUTPUT_DIRECTORY]

optional arguments:
  -h, --help            show this help message and exit
  -t TARGET, --target TARGET
                        the protein target (e.g. nsp13)
  -f FRAGMENTS [FRAGMENTS ...], --fragments FRAGMENTS [FRAGMENTS ...]
                        the list of fragments to merge
  -o OUTPUT_DIRECTORY, --output_directory OUTPUT_DIRECTORY
                        the directory to write the merge files to
```

**Example**
For example, to find all possible merges of fragments 34, 176 and 212, which are hits against nsp13, run the following: 

```
python scripts/run_database_query.py -t nsp13 -f x0034_0B x0176_0B x0212_0B -o data/example_folder
```
This will prompt you for your password to access the database. This will enumerate all possible pairs of fragments and for each pair, will query the database to find the synthons of fragment B, and then find expansions of fragment A using each synthon. For each pair, a json file is created in the specified output directory, named according to the merge, e.g. `x0034_0B_x0176_0B.json` will include all merges found involving expansions of fragment 34 using synthons of fragment 176. These are saved in a dictionary, with the dictionary keys indicating the synthon used in the expansion, and the values representing the list of merges found (in SMILES format).

N.B. Connecting to the database in this way is not always very stable, so I typically run this in batches, but this script shows how it is done.

### Filtering the merges

To filter the merges, `filtering.py` is run for each merge pair (i.e. each json file created according to the above). The options are shown below:

```
usage: filtering.py [-h] [-f MERGE_FILE] [-m MERGE] [-a FRAGMENT_A]
                    [-b FRAGMENT_B] [-p PROTEIN_A] [-q PROTEIN_B]
                    [-o OUTPUT_DIRECTORY]

optional arguments:
  -h, --help            show this help message and exit
  -f MERGE_FILE, --merge_file MERGE_FILE
                        the json file containing the merges
  -m MERGE, --merge MERGE
                        the name of the merge, e.g. x0107_0A_x0434_0A
  -a FRAGMENT_A, --fragment_A FRAGMENT_A
                        fragment A mol file
  -b FRAGMENT_B, --fragment_B FRAGMENT_B
                        fragment B mol file
  -p PROTEIN_A, --protein_A PROTEIN_A
                        protein pdb file associated with fragment A
  -q PROTEIN_B, --protein_B PROTEIN_B
                        protein pdb file associated with fragment B
  -o OUTPUT_DIRECTORY, --output_directory OUTPUT_DIRECTORY
                        the directory to write the filtered files to
```

**Example**
For example, to filter all of the merges of fragments 34 and 212, use the following:

```
python scripts/filtering.py -f data/example_folder/x0034_0B_x0212_0B.json -m x0034_0B_x0212_0B -a nsp13/aligned/nsp13-x0034_0B/nsp13-x0034_0B.mol -b nsp13/aligned/nsp13-x0212_0B/nsp13-x0212_0B.mol -p nsp13/aligned/nsp13-x0034_0B/nsp13-x0034_0B_apo-desolv.pdb -q nsp13/aligned/nsp13-x0212_0B/nsp13-x0212_0B_apo-desolv.pdb -o data/results_folder
```
Within the results directory a sub-directory called is created tempfiles, where all the Fragmenstein files are temporarily written to, and a Fragmenstein folder, to which the required Fragmenstein files (including the minimized mol and pdb files and the json containing the results needed for filtering) are saved. 

Running the filtering script creates two json files, a `MERGE_filtered.json` and a `MERGE_failures.json`. The former contains a dictionary containing the successfully filtered compounds, their unique name (which is used to name the Fragmenstein files), SMILES and PLIF score. The latter file contains the SMILES that did not pass the filters, their unique name and which filter they failed at. 

## Tests

The tests here require updating.
