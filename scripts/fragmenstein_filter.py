"""Place and filter the smiles with Fragmenstein"""

import pyrosetta
import json
import os
import shutil
import tempfile

from rdkit import Chem
from fragmenstein import Victor

def get_dict(json_file):
    """
    Function opens the json file to load the dictionary

    :param json_file: json containing the dictionary
    :type json_file: .json
    :return: dictionary containing Fragmenstein info
    :rtype: nested dictionary
    """
    with open(json_file) as f:
        data = json.load(f)
    return data

def create_directories(output_directory):
    """
    Create subdirectories to store the temporary files and fragmenstein results.
    """
    # create tempfiles folder
    if os.path.exists(os.path.join(output_directory, 'tempfiles')) == True:
        pass
    else:
        if not os.path.exists(output_directory):
            os.mkdir(output_directory)
        os.mkdir(os.path.join(output_directory, 'tempfiles'))

    # create fragmenstein folder to save the important files we want to keep
    if os.path.exists(os.path.join(output_directory, 'fragmenstein')) == True:
        pass
    else:
        os.mkdir(os.path.join(output_directory, 'fragmenstein'))

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
    # create temporary directory to write files to
    temp_dir = tempfile.TemporaryDirectory(dir=f'{output_directory}/tempfiles/')

    # run Fragmenstein
    # initialise PyRosetta
    pyrosetta.init(extra_options='-no_optH false -mute all -ex1 -ex2 -ignore_unrecognized_res false -load_PDB_components false -ignore_waters false')

    # get the fragments from the files
    fragments_fnames = [fragmentA, fragmentB]
    hits = [Chem.MolFromMolFile(frag) for frag in fragments_fnames]

    # set the output directory
    Victor.work_path = temp_dir.name

    # set up Victor and place the smiles
    v = Victor(hits=hits, pdb_filename=protein, covalent_resi= '2B')
    v.place(smiles=smiles, long_name=name)  # 'long_name' used to name the files

    # only keep the files needed for filtering (move to fragmenstein folder)
    # fragmenstein saves files with hyphens instead of underscores
    name_ = name.replace('_', '-')

    # move the json file
    json_file = name_ + '.minimised.json'
    json = os.path.join(temp_dir.name, name_, json_file)
    json_new = os.path.join(output_directory, 'fragmenstein', json_file)
    shutil.move(json, json_new)

    # move the mol file
    mol_file = name_ + '.minimised.mol'
    mol = os.path.join(temp_dir.name, name_, mol_file)
    mol_new = os.path.join(output_directory, 'fragmenstein', mol_file)
    shutil.move(mol, mol_new)

    # move the pdb file
    pdb_file = name_ + '.holo_minimised.pdb'
    pdb = os.path.join(temp_dir.name, name_, pdb_file)
    pdb_new = os.path.join(output_directory, 'fragmenstein', pdb_file)
    shutil.move(pdb, pdb_new)

    # remove files from temporary directory
    temp_dir.cleanup()

    return json_new, mol_new, pdb_new

def fragmenstein_filter(json_file):
    """
    Function filters molecules for those where both fragments were considered in its placement,
    a negative ΔΔG and combined RMSD with the fragments of < 1.5A.

    :param json_file: json containing the dictionary
    :type json_file: .json

    :return: returns 'pass' or 'fail'
    :rtype: string
    """
    data = get_dict(json_file) # load the dictionary from the json file

    # retrieve the energy of the bound and unbound molecules and the comRMSD
    G_bound = data['Energy']['ligand_ref2015']['total_score']  # energy of bound molecule
    G_unbound = data['Energy']['unbound_ref2015']['total_score']  # energy of unbound molecule
    deltaG = G_bound - G_unbound  # calculate energy difference
    comRMSD = data['mRMSD']  # RMSD between two fragments and merge

    # get number of fragments used for placement of SMILES
    regarded = 0
    for rmsd in data['RMSDs']:
            if rmsd != None:
                regarded += 1

    # only keep molecules where both fragments used in placement
    if regarded == 2:
        if deltaG < 0:  # keep molecules with negative ΔΔG
            if comRMSD <= 1:
                result = 'pass'
            else:
                result = 'fail'
        else:
            result = 'fail'  # remove molecule with positive ΔΔG
    else:
        result = 'fail'

    return result
