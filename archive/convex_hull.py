"""Elaboratability analysis"""
import argparse
import json
import os
from typing import Tuple

import numpy as np
from numpy import array
import pygamer
from joblib import Parallel, delayed
from pymol import cmd
from rdkit import Chem
from rdkit.Chem import AllChem, Mol, rdFMCS
from rdkit.Geometry.rdGeometry import Point3D
from scipy.spatial import ConvexHull
from tqdm import tqdm
from utils.filter_utils import sanitize
from utils.utils import get_distance

REACTION_INFO = {
    "Amidation": {
        "rxn": "[#6:1](=[#8:2])-[#8].[#7;H3,H2,H1:3]>>[#6:1](=[#8:2])-[#7:3]",
        "reactant1": Chem.MolFromSmiles("C(=O)O"),
        "reactant2": Chem.MolFromSmiles("CN"),
        "recursive_smarts1": Chem.MolFromSmarts("[O&$(OC(=O))]"),
        "recursive_smarts2": Chem.MolFromSmarts("[#7;H3,H2,H1:3]"),
    },
    "Amide_schotten-baumann": {
        "rxn": "[#7;H2,H1:3].[#6:1](=[#8:2])-[#17]>>[#6:1](=[#8:2])-[#7:3]",
        "reactant1": Chem.MolFromSmiles("CN"),
        "reactant2": Chem.MolFromSmiles("C(=O)Cl"),
        "recursive_smarts1": Chem.MolFromSmarts("[#7;H3,H2,H1:3]"),
        "recursive_smarts2": Chem.MolFromSmarts("[Cl&$(ClC(=O))]"),
    },
    "Reductive_amination": {
        "rxn": "[#6:2](=[#8])(-[#6:1]).[#7;H3,H2,H1:3]>>[#6:2](-[#6:1])-[#7:3]",
        "reactant1": Chem.MolFromSmiles("CC(=O)O"),
        "reactant2": Chem.MolFromSmiles("CN"),
        "recursive_smarts1": Chem.MolFromSmarts("[O&$(O=C(C))]"),
        "recursive_smarts2": Chem.MolFromSmarts("[#7;H3,H2,H1:3]"),
    },
    "N-nucleophilic_aromatic_substitution": {
        "rxn": "[#6:3]-[#7;H3,H2,H1:2].[c:1]-[F,Cl,Br,I]>>[#6:3]-[#7:2]-[c:1]",
        "reactant1": Chem.MolFromSmiles("CN"),
        "reactant2": Chem.MolFromSmiles("Fc1ccccc1"),
        "recursive_smarts1": Chem.MolFromSmarts("[N&$([#7;H3,H2,H1:2]-[#6:3])]"),
        "recursive_smarts2": Chem.MolFromSmarts("[F,Br,I,Cl&$([F,Br,I,Cl]c)]"),
    },
    "Sp2-sp2_Suzuki_coupling": {
        "rxn": "[c:1]-[F,Cl,Br,I].[#6:2]-[B]>>[c:1]-[#6:2]",
        "reactant1": Chem.MolFromSmiles("Fc1ccccc1"),
        "reactant2": Chem.MolFromSmiles("CB"),
        "recursive_smarts1": Chem.MolFromSmarts("[F,Br,I,Cl&$([F,Br,I,Cl]c)]"),
        "recursive_smarts2": Chem.MolFromSmarts("[B&$(Bc)]"),
    },
}


def get_products(mol: Mol) -> list:
    """
    Get all possible products using the set of defined reactions
    :param mol: mol to undergo elaboration
    :type mol: rdkit.Chem.Mol
    :return: list of elaborated molecules (products)
    :rtype: list
    """
    for atom in mol.GetAtoms():
        atom.SetIntProp("original_idx", atom.GetIdx())

    for bond in mol.GetBonds():
        bond.SetIntProp("original_idx", bond.GetIdx())

    reactants = ["reactant1", "reactant2"]

    products = []

    for reaction in REACTION_INFO:
        rxn = AllChem.ReactionFromSmarts(REACTION_INFO[reaction]["rxn"])
        for reactant in reactants:
            p1s = [
                p[0] for p in rxn.RunReactants((REACTION_INFO[reaction][reactant], mol))
            ]
            p2s = [
                p[0] for p in rxn.RunReactants((mol, REACTION_INFO[reaction][reactant]))
            ]
            products.extend(p1s)
            products.extend(p2s)

    return products



def get_vector_atoms(product: Mol, mol: Mol) -> Tuple[Mol, int, int, Point3D, str]:
    """
    For each product, use the MCS to get the attachment point atom and the atom it is bonded to in the elaboration

    :param product: elaborated product
    :type product: rdkit.Chem.Mol
    :param mol: mol that was elaborated
    :type mol: rdkit.Chem.Mol
    :return: tuple of product, product_exit_vector_idx, original_exit_vector_idx, product_exit_vector_coords, transformation_type
    :rtype: tuple
    """
    # load product molecule
    product = Chem.Mol(product)
    sanitize(product)

    # calculate MCS
    mcs = Chem.MolFromSmarts(rdFMCS.FindMCS([product, mol]).smartsString)
    sanitize(mcs)

    # use MCS to get which atoms have been added to the molecule in the reaction
    original_atom_idxs = list(product.GetSubstructMatch(mcs))
    elab_idxs = [atom.GetIdx() for atom in product.GetAtoms() if atom.GetIdx() not in original_atom_idxs]

    product_exit_vector_idx = None  # index of attachment point atom (for product)
    # product_added_atom_idx = None  # index of atom added directly to exit vector atom (for product)

    # using the mol bonds, get the bond that is between original and new atom (exit vector)
    # use this to get the relevant atom indices
    for bond in product.GetBonds():
        bond_idx = bond.GetIdx()
        bgn_idx = product.GetBondWithIdx(bond_idx).GetBeginAtomIdx()
        end_idx = product.GetBondWithIdx(bond_idx).GetEndAtomIdx()
        if bgn_idx in original_atom_idxs and end_idx in elab_idxs:
            product_exit_vector_idx = bgn_idx
            # product_added_atom_idx = end_idx

        elif bgn_idx in elab_idxs and end_idx in original_atom_idxs:
            product_exit_vector_idx = end_idx
            # product_added_atom_idx = bgn_idx

    # get the equivalent atom index in the original molecule (can be different to that in the product)
    original_exit_vector_idx = None
    product_exit_vector_coords = product.GetConformer().GetAtomPosition(product_exit_vector_idx)  # prod exit coords
    for atom in mol.GetAtoms():
        og_coords = list(mol.GetConformer().GetAtomPosition(atom.GetIdx()))
        if og_coords == list(product_exit_vector_coords):
            original_exit_vector_idx = atom.GetIdx()

    # get the transformation type (whether an original heavy atom has been replaced or new atoms attached)
    transformation_type = None
    if mcs.GetNumAtoms() < mol.GetNumAtoms():
        transformation_type = 'replacement'

    elif mcs.GetNumAtoms() == mol.GetNumAtoms():
        transformation_type = 'expansion'

    return product, product_exit_vector_idx, original_exit_vector_idx, product_exit_vector_coords, transformation_type


def get_attached_atom_coords_from_replacement(mol: Mol, product: Mol, product_exit_vector_idx: int, original_exit_vector_idx: int) -> np.array:
    """
    Get the coordinates of the atom attached to the exit vector atom (so the two sets of coordinates can be used as
    vector). The coordinates of the new attached atom will be completely different but the other attached atoms will
    be the same - so check the coordinates and get the one that doesn't match.

    :param mol: molecule to be elaborated
    :type mol: rdkit.Chem.Mol
    :param product: elaborated product
    :type product: rdkit.Chem.Mol
    :param product_exit_vector_idx: index of exit vector atom in the product
    :type product_exit_vector_idx: int
    :param original_exit_vector_idx: index of exit vector atom in the original molecule
    :type original_exit_vector_idx: int
    :return: numpy array of the exit vector coordinates
    :rtype: numpy.array
    """
    original_vector_attached_atoms = []

    for bond in mol.GetBonds():
        bond_idx = bond.GetIdx()
        bgn_idx = mol.GetBondWithIdx(bond_idx).GetBeginAtomIdx()
        end_idx = mol.GetBondWithIdx(bond_idx).GetEndAtomIdx()
        if bgn_idx == original_exit_vector_idx:
            if end_idx not in original_vector_attached_atoms:
                original_vector_attached_atoms.append(end_idx)
        if end_idx == original_exit_vector_idx:
            if bgn_idx not in original_vector_attached_atoms:
                original_vector_attached_atoms.append(bgn_idx)

    product_vector_attached_atoms = []

    for bond in product.GetBonds():
        bond_idx = bond.GetIdx()
        bgn_idx = product.GetBondWithIdx(bond_idx).GetBeginAtomIdx()
        end_idx = product.GetBondWithIdx(bond_idx).GetEndAtomIdx()
        if bgn_idx == product_exit_vector_idx:
            if end_idx not in product_vector_attached_atoms:
                product_vector_attached_atoms.append(end_idx)
        if end_idx == product_exit_vector_idx:
            if bgn_idx not in product_vector_attached_atoms:
                product_vector_attached_atoms.append(bgn_idx)

    original_vector_attached_coords = [list(mol.GetConformer().GetAtomPosition(idx)) for idx in original_vector_attached_atoms]
    product_vector_attached_coords = [list(product.GetConformer().GetAtomPosition(idx)) for idx in product_vector_attached_atoms]

    attached_coords = None
    for og_coords in original_vector_attached_coords:
        if og_coords not in product_vector_attached_coords:
            attached_coords = np.array(og_coords)

    return [attached_coords]


def get_attached_atom_coords_from_expansion(mol: Mol, mol_file: str, original_exit_vector_idx: int) -> list:
    """
    :param mol: molecule to be elaborated
    :type mol: rdkit.Chem.Mol
    :param mol_file: file path for molecule to be elaborated
    :type mol_file: str
    :param original_exit_vector_idx: index of exit vector atom in the original molecule
    :type original_exit_vector_idx: int
    :return: numpy array of the exit vector coordinates
    :rtype: numpy.array
    """
    mol_with_hs = Chem.MolFromMolFile(mol_file, removeHs=False)
    original_mol_atoms = [atom.GetIdx() for atom in mol.GetAtoms()]
    h_atoms = [atom.GetIdx() for atom in mol_with_hs.GetAtoms() if atom.GetIdx() not in original_mol_atoms]

    original_added_atom_idx = []
    for bond in mol_with_hs.GetBonds():
        bond_idx = bond.GetIdx()
        bgn_idx = mol_with_hs.GetBondWithIdx(bond_idx).GetBeginAtomIdx()
        end_idx = mol_with_hs.GetBondWithIdx(bond_idx).GetEndAtomIdx()
        if bgn_idx == original_exit_vector_idx and end_idx in h_atoms:
            original_added_atom_idx.append(end_idx)
        if end_idx == original_exit_vector_idx and bgn_idx in h_atoms:
            original_added_atom_idx.append(bgn_idx)

    all_added_atom_coords = [np.array(mol_with_hs.GetConformer().GetAtomPosition(i)) for i in original_added_atom_idx]
    return all_added_atom_coords


def get_vector_coords(mol, mol_file):
    """
    Get coordinates of the exit vector atoms (and the atoms in the elaboration they are attached to)
    """
    vectors = {}

    products = get_products(mol)
    for p in products:
        # get the indices of the exit vectors and the added atoms (for both the product and the original molecule)
        product, product_exit_vector_idx, original_exit_vector_idx, product_exit_vector_coords, transformation_type = get_vector_atoms(p, mol)

        if product_exit_vector_idx is not None and original_exit_vector_idx is not None:
            # check the exit vector is not already in the data dictionary
            if original_exit_vector_idx not in vectors.keys():
                try:
                    if transformation_type == 'replacement':
                        original_added_atom_coords = get_attached_atom_coords_from_replacement(mol, product, product_exit_vector_idx, original_exit_vector_idx)
                    if transformation_type == 'expansion':
                        original_added_atom_coords = get_attached_atom_coords_from_expansion(mol, mol_file, original_exit_vector_idx)

                    if len(original_added_atom_coords) > 0:
                        inner_d = {
                            "exit_coords": np.array(product_exit_vector_coords),
                            "added_atom_coords": original_added_atom_coords,  # list of np array
                        }
                        vectors[original_exit_vector_idx] = inner_d
                except:
                    print('Could not get added atom coordinates')
        else:
            print('No exit vector found')

    return vectors


def calc_angle_between_points(attach_coords, atom1_coords, atom2_coords):
    """
    Calculate the angle between three atoms (two elaboration vectors pointed at two different protein atoms) using
    their coordinates.
    """
    a = np.array(atom1_coords)
    b = np.array(attach_coords)
    c = np.array(atom2_coords)

    ba = a - b
    bc = c - b

    cos_angle = np.dot(ba, bc) / (np.linalg.norm(ba) * np.linalg.norm(bc))
    angle = np.arccos(cos_angle)
    return np.degrees(angle)


def remove_hydrogens(prot_file):
    """
    Remove hydrogens from pdb file before computing mesh
    """
    new_file = prot_file.replace('.pdb', '_no_hs.pdb')
    if not os.path.exists(new_file):
        cmd.reinitialize()
        cmd.load(prot_file, 'prot')
        cmd.select('elem H', 'hydrogens')
        cmd.remove('hydrogens')
        cmd.save(new_file)
    return new_file


def cone_analysis(
    attach_coords: np.array, all_added_atom_coords: list, protein_file: str, angle_limit=30, dist_lim=10
):
    hull_volumes = []
    prot_file_no_hs = remove_hydrogens(protein_file)
    mesh = pygamer.readPDB_molsurf(prot_file_no_hs)

    for added_atom_coords in all_added_atom_coords:
        vertex_coords = mesh.to_ndarray()[0]
        vertex_indices = []

        for i, vertex in enumerate(vertex_coords):
            angle = calc_angle_between_points(attach_coords, added_atom_coords, vertex)
            if angle <= angle_limit:
                vertex_indices.append(i)

        select_vertex_coords = vertex_coords[vertex_indices]

        if len(select_vertex_coords) == 0:
            print('No vertex coordinates within cone')
            return None

        dists = []
        for vertex_coord in select_vertex_coords:
            dist = get_distance(attach_coords, vertex_coord)
            dists.append(dist)

        dists.sort()
        if dists[0] > dist_lim:
            print('No vertex within', dist_lim)
            return None

        print('attach', attach_coords)
        print('added', added_atom_coords)
        attach_coords_ = attach_coords.reshape(1,3)
        added_atom_coords = added_atom_coords.reshape(1,3)

        try:
            all_coords = np.concatenate((attach_coords_, added_atom_coords, select_vertex_coords), axis=0)
            hull = ConvexHull(points=all_coords)
            volume = hull.volume
        except:
            volume = None

        hull_volumes.append(volume)

    return hull_volumes


def run_elab_analysis(i, mol_file, protein_file, output_dir, write_json=True):
    print(i)
    elab_file = protein_file.replace('.holo_minimised_nolig.pdb', '_cone_elab_hull.json')

    #TODO CHANGE BACK
    if not os.path.exists(elab_file) or os.path.exists(elab_file):

        try:
            mol = Chem.MolFromMolFile(mol_file)
            vectors = get_vector_coords(mol, mol_file)

            if len(vectors) > 0:
                for exit_vector in vectors:
                    # exit_coords = np.array(vectors[exit_vector]['exit_coords'])
                    # added_atom_coords = np.array(vectors[exit_vector]['added_atom_coords'])
                    exit_coords = vectors[exit_vector]["exit_coords"]
                    added_atom_coords = vectors[exit_vector]["added_atom_coords"]

                    hull_volumes = cone_analysis(
                        exit_coords, added_atom_coords, protein_file
                    )
                    vectors[exit_vector]["exit_coords"] = list(exit_coords)
                    vectors[exit_vector]["added_atom_coords"] = [list(coords) for coords in added_atom_coords]
                    vectors[exit_vector]["hull_volumes"] = hull_volumes

                if not write_json:
                    return vectors

                else:
                    with open(elab_file, "w") as f:
                        json.dump(vectors, f)

            else:
                if not write_json:
                    return None

                else:
                    with open(elab_file, "w") as f:
                        json.dump({}, f)
            return vectors
        except Exception as e:
            print(e, mol_file)
            with open(elab_file, "w") as f:
                json.dump({}, f)


def run_all_elab(dir):
    pairs = [i for i in os.listdir(dir) if os.path.isdir(os.path.join(dir, i))]
    pair_dirs = [os.path.join(dir, i) for i in pairs]

    # all_merge_dirs = []
    all_mol_files = []
    all_protein_files = []

    for pair_dir in pair_dirs:
        merges = [i for i in os.listdir(pair_dir) if os.path.isdir(os.path.join(pair_dir, i))]
        merge_dirs = [os.path.join(pair_dir, i) for i in merges]

        for merge, merge_dir in zip(merges, merge_dirs):
            mol_file = os.path.join(merge_dir, f"{merge}.minimised.mol")
            protein_file = os.path.join(merge_dir, f"{merge}.holo_minimised_nolig.pdb")

            if os.path.exists(protein_file) and os.path.exists(mol_file):
                all_mol_files.append(mol_file)
                all_protein_files.append(protein_file)

    print(len(all_mol_files), len(all_protein_files))
    idxs = [i for i in range(len(all_mol_files))]
    _ = Parallel(n_jobs=os.cpu_count(), backend="multiprocessing")(
        delayed(run_elab_analysis)(i, mol_file, protein_file)
        for i, mol_file, protein_file in tqdm(zip(idxs, all_mol_files, all_protein_files))
    )

#
# def main():
#     parser = argparse.ArgumentParser()
#     parser.add_argument("-d", "--dir", help="directory")
#     args = parser.parse_args()
#     run_all_elab(args.dir)
#
#
# if __name__ == "__main__":
#     main()
