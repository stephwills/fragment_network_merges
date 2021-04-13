import os
from rdkit import Chem
from fragmenstein import Victor
import pyrosetta

def place_with_fragmenstein(mergeNames, smiles, fragmentAs, fragmentBs, protein, output_directory):
    summaries = {}
    pyrosetta.init(extra_options='-no_optH false -mute all -ex1 -ex2 -ignore_unrecognized_res false -load_PDB_components false -ignore_waters false')  # initialise pyrosetta
    for i, (merge, smi, fragmentA, fragmentB) in enumerate(zip(mergeNames, smiles, fragmentAs, fragmentBs)):
        output_subdirectory = merge + '_' + str(i)
        fragments_fnames = [fragmentA, fragmentB]
        hits = [Chem.MolFromMolFile(frag) for frag in fragments_fnames]
        Victor.work_path = output_directory
        v = Victor(hits=hits, # list of rdkit molecules (your fragments)
                   pdb_filename=protein, # file name of apo protein
                   covalent_resi= '1A')
        v.place(smiles=smi,
                long_name=output_subdirectory,
                )
        v.make_pse()
        print("Placed merge", i)
        summaries[i]  = v.summarise()
        print(v.summarise())

    return summaries
