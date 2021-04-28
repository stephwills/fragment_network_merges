"""Place and filter the smiles with Fragmenstein"""

import pyrosetta
import json
import tempfile
from rdkit import Chem
from rdkit.Chem import rdmolfiles
from fragmenstein import Victor

def place_smiles(name, smiles, fragmentA, fragmentB, protein, output_directory):
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
    temp_dir = tempfile.TemporaryDirectory() # create temporary directory to write the output files to

    pyrosetta.init(extra_options='-no_optH false -mute all -ex1 -ex2 -ignore_unrecognized_res false -load_PDB_components false -ignore_waters false')  # initialise pyrosetta
    fragments_fnames = [fragmentA, fragmentB]  # filenames of the fragments
    hits = [Chem.MolFromMolFile(frag) for frag in fragments_fnames]  # gets the mols from the filenames
    Victor.work_path = temp_dir.name  # set the output directory
    v = Victor(hits=hits,  # list of rdkit molecules (fragments)
               pdb_filename=protein,  # file name of apo protein
               covalent_resi= '1A'
               )
    v.place(smiles=smiles,
            long_name=name,  # to name the files
            )

    # get the json file with the info we need from the temp directory and save it to the output folder
    # this is the only file we need for further filtering
    minimised_json = f'{temp_dir.name}/{name}/{name}.minimised.json'  # filepath
    minimised_json = minimised_json.replace('_', '-')
    json_object = get_dict(minimised_json)  # read in the json file
    new_filepath = f'{output_directory}/{name}.minimised.json'
    with open(new_filepath, 'w') as f:  # write to new file in permanent directory
        json.dump(json_object, f)

    # also read in the mol file and save it to the output folder (needed for interaction fp)
    minimised_mol_file = f'{temp_dir.name}/{name}/{name}.minimised.mol'  # filepath
    minimised_mol_file = minimised_mol_file.replace('_', '-')
    minimised_mol = rdmolfiles.MolFromMolFile(minimised_mol_file)  # read in the mol file
    new_mol_filepath = f'{output_directory}/{name}.minimised.mol'
    rdmolfiles.MolToMolFile(minimised_mol, new_mol_filepath)  # write to new file in permanent directory

    temp_dir.cleanup()  # remove files from temporary directory

    return new_filepath, new_mol_filepath

def get_dict(json_file):
    """
    Function opens the json file to load the dictionary

    :param json_file: json containing the dictionary
    :type json_file: .json
    """
    f = open(json_file)
    data = json.load(f)
    f.close()
    return data

def fragmenstein_filter(json_file):
    """
    Function filters molecules for those where both fragments were considered in its placement,
    a negative ΔΔG and combined RMSD with the fragments of < 1.5A.

    :param json_file: json containing the dictionary
    :type json_file: .json

    :return: returns 'pass' or 'fail'
    :rtype: string
    """
    data = get_dict(json_file) # load the dictionary from the json
    # retrieve the energy of the bound and unbound molecules and the comRMSD
    G_bound = data['Energy']['ligand_ref2015']['total_score']  # energy of bound molecule
    G_unbound = data['Energy']['unbound_ref2015']['total_score']  # energy of unbound molecule
    comRMSD = data['mRMSD']  # RMSD between two fragments and merge
    regarded = 0
    # get number of molecules regarded during placement
    for rmsd in data['RMSDs']:
            if rmsd != None:
                regarded += 1
    deltaG = G_bound - G_unbound  # calculate energy difference
    if regarded == 2:
        if deltaG < 0:  # keep molecules with negative ΔΔG
            if comRMSD <= 1:
                result = 'pass'
            else:
                result = 'fail'
    else:
        result = 'fail'
    return result
