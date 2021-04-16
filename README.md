# fragment_network_merges

Uses the fragment network to generate compounds that combine parts of two fragments. 
Merges are generated using 'synthon expansion' before undergoing numerous filtering steps to output a set of purchasable molecules that meet various criteria. 

Filtering steps involve filtering based on calculated descriptors, 3D filtering using constrained embedding (in which the atoms from the fragments are constrained using the original 3D coordinates), placing the merges in the protein and energy minimisation with Fragmenstein, and filtering based on the calculation of interaction fingerprints.

The notebooks demonstrate how compounds are generated and filtered using a test set of 12 fragment pairs.  

Currently the scripts are not suitable for their purpose and need refactoring. Many of the scripts have been written to process all SMILES rather than one single SMILES at a time. 
