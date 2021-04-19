"""Place the smiles with Fragmenstein"""
import pyrosetta
from rdkit import Chem
from fragmenstein import Victor

def place_with_fragmenstein(name, smiles, fragmentA, fragmentB, protein, output_directory):
    """
    Function places the merge in the protein using Fragmenstein. Minimisation is performed
    using pyrosetta. Creates an output folder containing .json file with metrics that
    can be used for filtering.

    :param name: name of the merge (to name the output subfolder/files)
    :type name: string
    :param smiles: smiles of the merge
    :type smiles: string
    :param fragmentA: fragment A file
    :type fragmentA: mol file
    :param fragmentB: fragment B file
    :type fragmentB: mol file
    :param protein: protein file
    :type protein: pdb file
    :param output_directory: filepath of output directory
    :type output_directory: filepath string
    """
    pyrosetta.init(extra_options='-no_optH false -mute all -ex1 -ex2 -ignore_unrecognized_res false -load_PDB_components false -ignore_waters false')  # initialise pyrosetta
    fragments_fnames = [fragmentA, fragmentB]  # filenames of the fragments
    hits = [Chem.MolFromMolFile(frag) for frag in fragments_fnames]  # gets the mols from the filenames
    Victor.work_path = output_directory  # set the output directory
    v = Victor(hits=hits,  # list of rdkit molecules (fragments)
               pdb_filename=protein,  # file name of apo protein
               covalent_resi= '1A'
               )
    v.place(smiles=smi,
            long_name=name,  # to name the files
            )
    v.make_pse()  # creates a pse file for pymol
