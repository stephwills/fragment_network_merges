This describes how to run the `overlap-filter` job from the `candidate generation` category in the `xchem` collection.

## What the job does

TODO - describe 

## Implementation details

* Job implementation: [overlap_filter.py](/filter/overlap_filter.py)
* Job definition: `jobs.overlap-filter` in [fragnet-merges.yaml](/data-manger/fragnet-merges.yaml)

## How to run the job

### Inputs

* **Compounds to check**: The molecules to test.
* **1st PDB file**: The first protein to check against.
* **2nd PDB file**: The second protein to check against.


### Options

* **Output file name**: Name for the output SD-file.
* **Distance filter for clashes**: Distance cutoff to protein in Angstroms .

### Outputs

A SD-file is output containing 3D structures that did not clash with the protein.
