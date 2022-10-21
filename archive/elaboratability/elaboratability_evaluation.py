"""Evaluate the elaboratability of compounds"""

from collections import OrderedDict

import numpy as np
import os
from joblib import Parallel, delayed
from pymol import cmd
from rdkit import Chem
from rdkit.Chem import (AllChem, rdFMCS, rdmolops, rdMolTransforms,
                        rdShapeHelpers)
from scipy import stats
from tqdm import tqdm
from utils.filter_utils import add_coordinates
from utils.utils import get_distance

REACTIONS = {
    "Amidation": "[#6:1](=[#8:2])-[#8].[#7;H3,H2,H1:3]>>[#6:1](=[#8:2])-[#7:3]",
    "Amide_schotten-baumann": "[#7;H2,H1:3].[#6:1](=[#8:2])-[#17]>>[#6:1](=[#8:2])-[#7:3]",
    "Reductive_amination": "[#6:2](=[#8])(-[#6:1]).[#7;H3,H2,H1:3]>>[#6:2](-[#6:1])-[#7:3]",
    "N-nucleophilic_aromatic_substitution": "[#6:3]-[#7;H3,H2,H1:2].[c:1]-[F,Cl,Br,I]>>[#6:3]-[#7:2]-[c:1]",
    "Sp2-sp2_Suzuki_coupling": "[c:1]-[F,Cl,Br,I].[#6:2]-[B]>>[c:1]-[#6:2]",
}

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

SMILES = [
    "C(=O)O",
    "CN",
    "C(=O)Cl",
    "CC(=O)O",
    "Fc1ccccc1",
    "CB",
]  # fragments cover possible reactants for reactions


FRAGMENTS = [
    "F",
    "O",
    "S",
    "C=O",
    "C=C",
    "C#C",
    "C",
    "C1CC1",
    "C1CCC1",
    "C1CCCC1",
    "C1C=CC=C1",
    "C1CCCCC1",
    "C(c1ccccc1)",  # check; doesn't work directly with aromatic ring
    "C1CCCCCC1",
    "N",
    "N(=O)O",
    "P",
    "P(=O)(O)(O)O",
    "S(=O)(=O)(O)O",
    "C(=O)O",
    "C(=O)N",
    "N",
    "S(=O)=O",
    "S(=O)(=O)N",
]


def run_reactions(mol, reactants):
    """
    Run reactions for one filtered compound against all reaction types and reactants.
    Returns list of products and the substruct match with the original molecule (so you can see what's been retained)
    """
    products = []
    chosen_reactants = []
    patts = []

    for reaction in REACTIONS:
        rxn = AllChem.ReactionFromSmarts(REACTIONS[reaction])
        ps_patt = REACTIONS[reaction].split(">>")[1]
        for reactant in reactants:
            ps1 = rxn.RunReactants((mol, reactant))
            ps2 = rxn.RunReactants((reactant, mol))

            for ps in [ps1, ps2]:
                if ps:
                    for p in ps:
                        products.append(p[0])
                        patts.append(
                            p[0].GetSubstructMatch(Chem.MolFromSmarts(ps_patt))
                        )
                        reactants.append(reactant)

    return products, patts, chosen_reactants


def replace_hydrogens(mol, substructs):
    """
    Replace all hydrogens with selected functional group
    """
    products = []
    chosen_substructs = []
    attch_idxs = []

    for substruct in substructs:
        substruct = Chem.MolFromSmiles(substruct)
        mol = Chem.AddHs(
            Chem.RWMol(Chem.MolFromSmiles(Chem.MolToSmiles(mol))), addCoords=True
        )
        atoms = [atom.GetIdx() for atom in mol.GetAtoms()]
        repl = rdmolops.ReplaceSubstructs(mol, Chem.MolFromSmarts("[H]"), substruct)

        for rep in repl:
            mcs = rdFMCS.FindMCS([mol, rep]).smartsString
            mcs_mol = Chem.MolFromSmarts(mcs)
            attch_idx = [i for i in atoms if i not in mol.GetSubstructMatch(mcs_mol)]
            attch_idxs.append(attch_idx)
            chosen_substructs.append(substruct)
            products.append(rep)

    return products, attch_idxs, chosen_substructs


def embed_elab(mol, elab, addHs=False):
    """
    Perform constrained embedding on the elaborated molecule
    """
    _ = AllChem.Compute2DCoords(elab)

    # calculates maximum common substructure (MCS) between the fragment and the elaborated compound
    mcs = rdFMCS.FindMCS([elab, mol])
    mcs_mol = Chem.MolFromSmarts(mcs.smartsString)

    # add coordinates to the MCS
    matches = mol.GetSubstructMatch(mcs_mol)
    core = add_coordinates(mol, mcs_mol, matches)

    # perform constrained embedding (may fail for some molecules)
    nm = Chem.Mol(elab)
    if addHs:
        try:
            nm = Chem.AddHs(nm)  # add hydrogens for more accurate embedding
        except:
            print("Molecule already has hydrogens")
    embedded = AllChem.ConstrainedEmbed(nm, core)
    AllChem.UFFOptimizeMolecule(embedded)
    embedded = Chem.RemoveHs(embedded)

    return embedded


def run_reaction_eval(mol, prot, reactants=FRAGMENTS, n_confs=100):
    """
    Evaluate elaboratability using reactions for a molecule
    """
    products, patts, chosen_reactants = run_reactions(mol, reactants)
    clashes = []
    for product in products:
        clash = avg_clash(mol, product, prot, n_confs)
        clashes.append(clash)

    return products, patts, chosen_reactants, clashes


def run_hydrogens_eval(mol, prot, substructs=FRAGMENTS, n_confs=100):
    """
    Evaluate elaboratability using hydrogen replacement for a molecule
    """
    products, attch_idxs, chosen_substructs = replace_hydrogens(mol, substructs)
    clashes = []
    for product in tqdm(products):
        clash = avg_clash(mol, product, prot, n_confs)
        clashes.append(clash)

    return products, attch_idxs, chosen_substructs, clashes


def avg_clash_rotate(mol, elab, prot, rotat_int, returnRotated=False):
    """
    For a given elaboration, rotate the torsional angle joining the functional group and calculate the average clash
    with the protein
    """
    product = Chem.Mol(elab)
    mcs = Chem.MolFromSmarts(rdFMCS.FindMCS([product, mol]).smartsString)
    matches = product.GetSubstructMatch(mcs)
    non_matches = [i.GetIdx() for i in product.GetAtoms() if i.GetIdx() not in matches]

    connecting_neighbours = None
    atom3 = None

    for non_match in non_matches:
        for at in product.GetAtomWithIdx(non_match).GetNeighbors():
            if at.GetIdx() not in non_matches:
                atom3 = non_match  # connecting atom
                connecting_neighbours = [
                    i.GetIdx() for i in product.GetAtomWithIdx(non_match).GetNeighbors()
                ]

    atom1 = [i for i in matches if i not in connecting_neighbours][
        0
    ]  # atom from original molecule (not neighbour)
    atom2 = [i for i in connecting_neighbours if i in matches][
        0
    ]  # atom from original molecule (neighbour)
    atom4 = [i for i in connecting_neighbours if i in non_matches][
        0
    ]  # atom from elab (not the connector)

    embedded_product = embed_elab(mol, product)

    rotated = []
    dists = []
    for interval in range(0, 360, rotat_int):
        new_elab = Chem.Mol(embedded_product)
        rdMolTransforms.SetDihedralDeg(
            new_elab.GetConformer(), atom1, atom2, atom3, atom4, interval
        )
        rotated.append(new_elab)
        dist = 1 - rdShapeHelpers.ShapeProtrudeDist(new_elab, prot)
        dists.append(dist)

    mean = np.mean(dists)

    if returnRotated:
        return mean, rotated

    else:
        return mean


def run_rotat_eval(mol, prot, substructs, rotat_int):
    """
    Add certain groups and generate conformations by instead rotating the torsional angle connecting the functional
    group
    """
    products, attch_idxs, chosen_substructs = replace_hydrogens(mol, substructs)

    clashes = []
    for product in products:
        clash = avg_clash_rotate(mol, product, prot, rotat_int)
        clashes.append(clash)

    return products, attch_idxs, chosen_substructs, clashes


# may be useful later
def _avg_clash_rotate(elab, prot, atom1, atom2, atom3, atom4, rot_int=45):
    """
    Rotate torsion angle and calculate avg clash distance of each conformation
    """
    dists = []
    for interval in range(0, 360, rot_int):
        new_elab = Chem.Mol(elab)
        rdMolTransforms.SetDihedralDeg(
            new_elab.GetConformer(), atom1, atom2, atom3, atom4, interval
        )
        dist = 1 - rdShapeHelpers.ShapeProtrudeDist(new_elab, prot)
        dists.append(dist)

    mean = np.mean(dists)
    return mean


def _get_rotated(mol, elab, rotat_int):
    product = Chem.Mol(elab)
    mcs = Chem.MolFromSmarts(rdFMCS.FindMCS([product, mol]).smartsString)
    matches = product.GetSubstructMatch(mcs)
    non_matches = [i.GetIdx() for i in product.GetAtoms() if i.GetIdx() not in matches]

    connecting_neighbours = None
    atom3 = None

    for non_match in non_matches:
        for at in product.GetAtomWithIdx(non_match).GetNeighbors():
            if at.GetIdx() not in non_matches:
                atom3 = non_match  # connecting atom
                connecting_neighbours = [
                    i.GetIdx() for i in product.GetAtomWithIdx(non_match).GetNeighbors()
                ]
    # atom from original molecule (not neighbour)
    atom1 = [i for i in matches if i not in connecting_neighbours][0]
    # atom from original molecule (neighbour)
    atom2 = [i for i in connecting_neighbours if i in matches][0]
    # atom from elab (not the connector)
    atom4 = [i for i in connecting_neighbours if i in non_matches][0]

    embedded_product = embed_elab(mol, product)

    rotated = []
    for interval in range(0, 360, rotat_int):
        new_elab = Chem.Mol(embedded_product)
        rdMolTransforms.SetDihedralDeg(
            new_elab.GetConformer(), atom1, atom2, atom3, atom4, interval
        )
        rotated.append(new_elab)

    return rotated


def dist_to_closest_atoms(mol, atom_idx, prot, n_dists=5):
    """
    Evaluate the mean distance to the N number of closest atoms in the protein to a selected elaboration vector
    (denoted using the atom index).
    """
    conf = mol.GetConformer()
    pos = np.array(
        conf.GetAtomPosition(atom_idx)
    )  # get coordinates of the ligand exit vector atom

    dists = []
    prot_conf = prot.GetConformer()
    for i in range(prot_conf.GetNumAtoms()):
        prot_pos = np.array(
            prot_conf.GetAtomPosition(i)
        )  # get the positions of all protein atoms
        dist = get_distance(
            pos, prot_pos
        )  # calculate the distance away from the ligand atom
        dists.append(dist)

    # take the geometric mean of the closest distances
    dists.sort()
    closest_dists = dists[:n_dists]
    mean = stats.mstats.gmean(closest_dists)
    return mean

#
# def dist_to_closest_spread_atoms(mol, atom_idx, prot, n_dists=5, min_angle=10):
#     """
#     Evaluate the mean distance to the N number of closest atoms in the protein to a selected elaboration vector
#     (denoted using the atom index).
#     """
#
#     conf = mol.GetConformer()
#     pos = np.array(
#         conf.GetAtomPosition(atom_idx)
#     )  # get coordinates of the ligand exit vector atom
#
#     dists = []
#     prot_coords = []
#
#     prot_conf = prot.GetConformer()
#     for i in range(prot_conf.GetNumAtoms()):
#         prot_pos = np.array(
#             prot_conf.GetAtomPosition(i)
#         )  # get the positions of all protein atoms
#         prot_coords.append(prot_pos)
#         dist = get_distance(
#             pos, prot_pos
#         )  # calculate the distance away from the ligand atom
#         dists.append(dist)
#
#     # take the geometric mean of the closest distances
#     prot_coords = [coord for _,coord in sorted(zip(dists, prot_coords), key=lambda x: x[0])]
#     dists.sort()
#
#     # record the closest 'spread out' dists (and the coords of the protein atoms)
#     spread_dists = []
#     spread_coords = []
#
#     # keep iterating over the closest dists to get those that are not in the same direction
#     c = 0
#
#     for coord, dist in zip(prot_coords, dists):
#         if len(spread_dists) == 0:
#             spread_dists.append(dist)
#             spread_coords.append(coord)
#             c += 1
#
#         else:
#             if c < n_dists:
#                 angles = []
#                 for spread_dist, spread_coord in zip(spread_dists, spread_coords):
#                     new_coord = coord
#                     angle = calc_angle_between_points(pos, spread_coord, new_coord)
#                     angles.append(angle)
#
#                 # if the min angle against the already recorded atoms/dists is < the min angle, then add to the list
#                 # print(angles, dist)
#                 if min(angles) >= min_angle:
#                     spread_dists.append(dist)
#                     spread_coords.append(coord)
#                     c += 1
#
#     mean = stats.mstats.gmean(spread_dists)
#     return mean


def dist_to_closest_spread_atoms(mol, atom_idx, protein_file, n_dists=5, min_angle=10):
    """
    Evaluate the mean distance to the N number of closest atoms in the protein to a selected elaboration vector
    (denoted using the atom index).
    """

    conf = mol.GetConformer()
    pos = np.array(
        conf.GetAtomPosition(atom_idx)
    )  # get coordinates of the ligand exit vector atom

    dists = []
    prot_coords = []


    cmd.reinitialize()
    cmd.load(protein_file, 'protein')
    prot_coords = cmd.get_coords('all')
    for prot_pos in prot_coords:
        dist = get_distance(pos, prot_pos)
        dists.append(dist)

    # take the geometric mean of the closest distances
    prot_coords = [coord for _,coord in sorted(zip(dists, prot_coords), key=lambda x: x[0])]
    dists.sort()

    # record the closest 'spread out' dists (and the coords of the protein atoms)
    spread_dists = []
    spread_coords = []

    # keep iterating over the closest dists to get those that are not in the same direction
    c = 0

    for coord, dist in zip(prot_coords, dists):
        if len(spread_dists) == 0:
            spread_dists.append(dist)
            spread_coords.append(coord)
            c += 1

        else:
            if c < n_dists:
                angles = []
                for spread_dist, spread_coord in zip(spread_dists, spread_coords):
                    new_coord = coord
                    angle = calc_angle_between_points(pos, spread_coord, new_coord)
                    angles.append(angle)

                # if the min angle against the already recorded atoms/dists is < the min angle, then add to the list
                # print(angles, dist)
                if min(angles) >= min_angle:
                    spread_dists.append(dist)
                    spread_coords.append(coord)
                    c += 1

    mean = stats.mstats.gmean(spread_dists)
    return mean


def _dists_to_closest_atoms(mol_file, prot_file, spread=True, n_dists=5, min_angle=10):
    """
    Get all attachment indices and mean pocket distance for a given mol file and protein file
    """
    try:
        # repair_pdb = prot_file.replace('.holo_minimised_nolig.pdb', '_Repair.pdb')
        # if os.path.exists(repair_pdb):
        #     prot_file = repair_pdb
        mol = Chem.MolFromMolFile(mol_file)
        # prot = Chem.MolFromPDBFile(prot_file)
        attach_idxs = get_attachment_points(mol)
        dists = []
        for idx in attach_idxs:
            if spread:
                dists.append(dist_to_closest_spread_atoms(mol, idx, prot_file, n_dists, min_angle))
            # else:
            #     dists.append(dist_to_closest_atoms(mol, idx, prot, n_dists))
        return attach_idxs, dists
    except Exception as e:
        print(mol_file, e)


def calc_dists_to_closest_atoms(mol_files, prot_files, spread=True, n_dists=10, min_angle=10):
    """
    Return attachment index and pocket distance calcualtion for all mol and protein files in parallel
    """
    res = Parallel(n_jobs=os.cpu_count(), backend="multiprocessing")(
        delayed(_dists_to_closest_atoms)(mol_file, prot_file, spread, n_dists, min_angle)
        for mol_file, prot_file in tqdm(zip(mol_files, prot_files))
    )
    all_idxs = [r[0] for r in res]
    all_dists = [r[1] for r in res]
    return all_idxs, all_dists


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


def get_attachment_points(mol, reaction_dict=REACTION_INFO):
    """
    Returns the indices of attachment point atoms by checking for substructure match with the reactants for some
    specified reactions. The index of the exact attachment point atom is identified by using recursive SMARTS (these
    have been written out manually).
    """
    attach_idxs = set()
    patts = ["recursive_smarts1", "recursive_smarts2"]
    for reaction in reaction_dict:
        for patt in patts:
            matches = mol.GetSubstructMatches(reaction_dict[reaction][patt])
            if len(matches) > 0:
                # print(reaction_dict[reaction][patt], matches)
                for match in matches:
                    attach_idxs.add(match[0])
    return list(attach_idxs)


def get_attachment_counts(mols, reaction_dict=REACTION_INFO):
    """
    Count attachment points
    """
    all_attach_idxs = []
    for mol in mols:
        attach_idxs = get_attachment_points(mol)
        all_attach_idxs.append(attach_idxs)
    counts = [len(i) for i in all_attach_idxs]
    return all_attach_idxs, counts


### second half of pipeline

def add_group(mol, idx, mod, mod_pos=0):
    """
    Add functional group to molecule given the attachment index, the smiles of the group to modify it with
    and the relative position within the group of the atom to bond to the molecule (is 0 if the SMILES is ordered
    with this atom first).
    """
    # get attachment atom
    attach_atom = mol.GetAtomWithIdx(idx)

    # count the number of neighbours of the attachment index
    neighs = attach_atom.GetNeighbors()
    n_neighs = len(neighs)

    rwmol = Chem.RWMol(mol)

    # get the max index in the mol
    mol_idxs = [i.GetIdx() for i in mol.GetAtoms()]
    max_idx = max(mol_idxs)

    if n_neighs == 1:
        # get idx of the neighbor
        neigh_idx = neighs[0].GetIdx()
        rwmol.RemoveAtom(idx)
        bond_idx = neigh_idx
        mod_atom_idx = max_idx + mod_pos

    else:
        bond_idx = idx
        # get the atom within the functional group to bond to the molecule
        # the new atoms are numbered from the max mol atom, then use relative position of atom to get new idx
        mod_atom_idx = max_idx + 1 + mod_pos

    # print(mod_atom_idx)
    combo = Chem.CombineMols(rwmol, Chem.MolFromSmiles(mod))
    ed_mol = Chem.RWMol(combo)
    ed_mol.AddBond(bond_idx, mod_atom_idx, order=Chem.rdchem.BondType.SINGLE)
    return_mol = Chem.Mol(ed_mol)
    # return_mol = bad_sanitize(return_mol)

    return return_mol




def embed_multiple_elab(mol, elab, addHs=False, n_confs=100):
    _ = AllChem.Compute2DCoords(elab)

    # calculates maximum common substructure (MCS) between the fragment and the elaborated compound
    mcs = rdFMCS.FindMCS([elab, mol])
    mcs_mol = Chem.MolFromSmarts(mcs.smartsString)
    mcs_mol = bad_sanitize(mcs_mol)
    # add coordinates to the MCS
    matches = mol.GetSubstructMatch(mcs_mol)
    core = add_coordinates(mol, mcs_mol, matches)

    embedded_mols = []
    for n in range(n_confs):
        # perform constrained embedding (may fail for some molecules)
        nm = Chem.Mol(elab)
        if addHs:
            try:
                nm = Chem.AddHs(nm)  # add hydrogens for more accurate embedding
            except:
                print("Molecule already has hydrogens")
        embedded = AllChem.ConstrainedEmbed(nm, core)
        AllChem.UFFOptimizeMolecule(embedded)
        embedded = Chem.RemoveHs(embedded)
        embedded_mols.append(embedded)

    return embedded_mols


def avg_clash(mol, elab, prot, n_confs=100, addHs=False):
    """
    Perform multiple constrained embeddings and calculate the avg clash distance
    """
    dists = []
    embedded_mols = embed_multiple_elab(mol, elab, n_confs=n_confs, addHs=addHs)

    for embedded in embedded_mols:
        dist = 1 - rdShapeHelpers.ShapeProtrudeDist(embedded, prot)
        dists.append(dist)

    mean = np.mean(dists)
    return mean


def add_all_groups(mol, idxs, mods):
    """
    Get all modified mols given a set of attachment indixes and the smiles of the functional groups you want to add.
    """
    all_modified_mols = []
    for idx in idxs:
        modified_mols = []
        for mod in mods:
            mod_mol = add_group(mol, idx, mod)
            modified_mols.append(mod_mol)

        all_modified_mols.append(modified_mols)

    return all_modified_mols


FUNCTIONAL_GROUPS = [
    "C",
    "N",
    "C(=O)O",
    "C1CCCCC1"
]

def bad_sanitize(mol):
    try:
        Chem.SanitizeMol(mol)
    except:
        mol.UpdatePropertyCache(strict=False)
        Chem.SanitizeMol(
            mol,
            Chem.SanitizeFlags.SANITIZE_FINDRADICALS
            | Chem.SanitizeFlags.SANITIZE_KEKULIZE
            | Chem.SanitizeFlags.SANITIZE_SETAROMATICITY
            | Chem.SanitizeFlags.SANITIZE_SETCONJUGATION
            | Chem.SanitizeFlags.SANITIZE_SETHYBRIDIZATION
            | Chem.SanitizeFlags.SANITIZE_SYMMRINGS,
            catchErrors=True,
        )
    return mol

def secondary_pipeline(mol, prot, idxs, mods=FUNCTIONAL_GROUPS, addHs=False, n_confs=100):
    clash_dists = {}

    for idx in idxs:
        idx_dists = {}
        for mod in mods:
            mod_mol = add_group(mol, idx, mod)
            # embedded_mols = embed_multiple_elab(mol, mod_mol, addHs=addHs, n_confs=n_confs)
            mean_clash = avg_clash(mol, mod_mol, prot, n_confs=n_confs, addHs=addHs)
            idx_dists[mod] = mean_clash
        clash_dists[idx] = idx_dists

    return clash_dists








"""
Pipeline

For molecule - generate attachment idxs
Calc distances
Get best attach idxs/molecules

Input - Molecule and set of indices
For each index - add the functional groups

{idx: {group: elab}}

{idx: {group: mean clash,
       group: mean clash}
    }
    


"""
