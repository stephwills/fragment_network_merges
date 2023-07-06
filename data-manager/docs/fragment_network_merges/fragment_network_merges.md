This describes how to run the `fragment_network-merge` jobs from the `comp chem`
category in the `fragment_network` collection.

## What the job does

fragment_network_merges provides a pipeline for querying the Fragment Network to find fragment merges and 
filtering the results to find the most promising compound. First, 2D candidante are pulled from a neo4j database.
Then, we filter the merges down to a manageable number of the most promising compounds.

See the [fragment_network_merges documentation](https://github.com/stephwills/fragment_network_merges) for more information


## Implementation details

* Job implementation: [fragment_network-merge](/frag_merge.nf) #TODO: Not sure nf required for such a simple job
* Job definition: `jobs.fragment_network-merge` in [fragment_network.yaml](../fragment_network.yaml),

## How to run the job

### Inputs

* **Fragment molecules**: Two or more fragment molecules as molfile or SD-file (1).
* **PDB file for protein**: The protein(s) in which the generated molecules are minimised. Ideally, you should provide one protein for each molecule

Notes:
(1) Input fragments are specified either as separate molfiles or a single SD-file. Specify a relatively small number of
fragments (e.g. less than 20).

### Options

* **Output file name**: Name for the output SD-file.
* **Input field name containing the fragment ID**: Field in the input fragments that will be written to the outputs as
the `ref_mols` field.
* **Include SMILES in output using this field name**: If specified the SMILES of the generated molecule is written to
a field of this name in the output molecules.
* **Include PDB details in output using this field name**: If specified details of the PDB file are written to a field
of this name in the output molecules.
* **Use this value for the proteinFieldName**: Use this value for the proteinFieldName value. If not specified the PDB file name is used.



### Outputs

A SD-file is output containing 3D structures for the generated molecules. *Number of molecules to generate* molecules are
generated for each set of input fragments. Some input fragments may not generate any outputs.

The following fields are written to the SD-file:

* **IDX**: The index of the molecule (related to *Number of molecules to generate*).
* **DDG**: The difference in delta G between the bound and unbound minimised molecule.
* **RMSD** The mean RMSD between the fragments and the minimised molecule.
* **Value of "Includes SMILES in output using this field name"**: If this option is specified the SMILES of the generated
  molecule is written to this field
**Value of "Include PDB details in output using this field name"**: If this option is specified details of the PDB file 
are written to this field. The value written is the value of the "Use this value for the proteinFieldName" property if
defined, or the PDB file name if not.
* **ref_mols**: If the *Input field name containing the fragment ID* option is specified then the value of this field is
used as the ID of the input fragment and these IDs are written to this field as a comma separated list (3).
#TODO: Check which scores are used by Steph
* **sa_score**: The synthetic accessibility score. 1 means easy, 10 means hard.
* **SuCOS_Score**: SuCOS score (0=none, 1=perfect)
* **SuCOS_FeatureMap_Score**: SuCOS feature map score (0=none, 1=perfect)
* **SuCOS_Protrude_Score**: SuCOS protrude score (0=none, 1=perfect)
* **o3da_score**: Open3DAlign score (larger numbers indicate better alignment)
* **o3da_score_rel**: Open3DAlign relative score (0=none, 1=perfect)
* **o3da_align**: Open3DAlign alignment score (larger numbers indicate better alignment)

Notes:
(3) e.g. if the input fragments have their IDs as their title line (first line in the record) then you can specify 
`_Name` as the value for this field and those IDs will be written to the `ref_mols` field e.g. as `Mpro-x0072_0A,Mpro-x0104_0A`.
If the input is a SD-file then you can also specify the name of any regular SD-file field. If you do not specify a value
for this option then the `ref_mols` field is not written.
