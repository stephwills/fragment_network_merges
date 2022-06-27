"""Evaluate the elaboratability of compounds"""

import numpy as np
from rdkit import Chem
from rdkit.Chem import (AllChem, rdFMCS, rdmolops, rdMolTransforms,
                        rdShapeHelpers)
from utils.filter_utils import add_coordinates

REACTIONS = {
    "Amidation": "[#6:1](=[#8:2])-[#8].[#7;H3,H2,H1:3]>>[#6:1](=[#8:2])-[#7:3]",
    "Amide_schotten-baumann": "[#7;H2,H1:3].[#6:1](=[#8:2])-[#17]>>[#6:1](=[#8:2])-[#7:3]",
    "Reductive_amination": "[#6:2](=[#8])(-[#6:1]).[#7;H3,H2,H1:3]>>[#6:2](-[#6:1])-[#7:3]",
    "N-nucleophilic_aromatic_substitution": "[#6:3]-[#7;H3,H2,H1:2].[c:1]-[F,Cl,Br,I]>>[#6:3]-[#7:2]-[c:1]",
    "Sp2-sp2_Suzuki_coupling": "[c:1]-[F,Cl,Br,I].[#6:2]-[B]>>[c:1]-[#6:2]",
}

SMILES = ["C(=O)O", "CN", "C(=O)Cl", "CC(=O)O", "Fc1ccccc1", "CB"]


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


def embed_elab(mol, elab):
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
    try:
        nm = Chem.AddHs(nm)  # add hydrogens for more accurate embedding
    except:
        print("Molecule already has hydrogens")
    embedded = AllChem.ConstrainedEmbed(nm, core)
    AllChem.UFFOptimizeMolecule(embedded)
    embedded = Chem.RemoveHs(embedded)

    return embedded


def avg_clash(mol, elab, prot, n_confs=100):
    """
    Perform multiple constrained embeddings and calculate the avg clash distance
    """
    dists = []
    for n in range(n_confs):
        embedded = embed_elab(mol, elab)
        dist = 1 - rdShapeHelpers.ShapeProtrudeDist(embedded, prot)
        dists.append(dist)

    mean = np.mean(dists)
    return mean


def run_reaction_eval(mol, prot, reactants=SMILES, n_confs=100):
    """
    Evaluate elaboratability using reactions for a molecule
    """
    products, patts, chosen_reactants = run_reactions(mol, reactants)
    clashes = []
    for product in products:
        clash = avg_clash(mol, product, prot, n_confs)
        clashes.append(clash)

    return products, patts, chosen_reactants, clashes


def run_hydrogens_eval(mol, prot, substructs=SMILES, n_confs=100):
    """
    Evaluate elaboratability using hydrogen replacement for a molecule
    """
    products, attch_idxs, chosen_substructs = replace_hydrogens(mol, substructs)
    clashes = []
    for product in products:
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
                connecting_neighbours = [i.GetIdx() for i in product.GetAtomWithIdx(non_match).GetNeighbors()]

    atom1 = [i for i in matches if i not in connecting_neighbours][0]  # atom from original molecule (not neighbour)
    atom2 = [i for i in connecting_neighbours if i in matches][0]  # atom from original molecule (neighbour)
    atom4 = [i for i in connecting_neighbours if i in non_matches][0]  # atom from elab (not the connector)

    embedded_product = embed_elab(mol, product)

    rotated = []
    dists = []
    for interval in range(0, 360, rotat_int):
        new_elab = Chem.Mol(embedded_product)
        rdMolTransforms.SetDihedralDeg(new_elab.GetConformer(), atom1, atom2, atom3, atom4, interval)
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
                connecting_neighbours = [i.GetIdx() for i in product.GetAtomWithIdx(non_match).GetNeighbors()]

    atom1 = [i for i in matches if i not in connecting_neighbours][0]  # atom from original molecule (not neighbour)
    atom2 = [i for i in connecting_neighbours if i in matches][0]  # atom from original molecule (neighbour)
    atom4 = [i for i in connecting_neighbours if i in non_matches][0]  # atom from elab (not the connector)

    embedded_product = embed_elab(mol, product)

    rotated = []
    for interval in range(0, 360, rotat_int):
        new_elab = Chem.Mol(embedded_product)
        rdMolTransforms.SetDihedralDeg(new_elab.GetConformer(), atom1, atom2, atom3, atom4, interval)
        rotated.append(new_elab)

    return rotated
