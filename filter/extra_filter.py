"""Extra filters"""

from rdkit import Chem
from rdkit.Chem.Scaffolds import MurckoScaffold
from rdkit.Chem.rdMolDescriptors import CalcNumHeavyAtoms, CalcNumRings
from rdkit.Chem import rdmolops, AllChem


# code from https://iwatobipen.wordpress.com/2020/01/23/cut-molecule-to-ring-and-linker-with-rdkit-rdkit-chemoinformatics-memo/


def is_in_samering(idx1, idx2, bond_rings):
    for bond_ring in bond_rings:
        if idx1 in bond_ring and idx2 in bond_ring:
            return True
    return False


def getLinkerbond(mol, useScaffold=True):
    res = []
    for atom in mol.GetAtoms():
        atom.SetIntProp("orig_idx", atom.GetIdx())
    for bond in mol.GetBonds():
        bond.SetIntProp("orig_idx", bond.GetIdx())

    if useScaffold:
        mol = MurckoScaffold.GetScaffoldForMol(mol)

    ring_info = mol.GetRingInfo()
    bond_rings = ring_info.BondRings()  # gives tuples with bonds that are in the rings
    ring_bonds = set()
    for ring_bond_idxs in bond_rings:  # add all the bonds that are in rings to a set
        for idx in ring_bond_idxs:
            ring_bonds.add(idx)
    all_bonds_idx = [bond.GetIdx() for bond in mol.GetBonds()]  # get all bond idxs in list
    none_ring_bonds = set(all_bonds_idx) - ring_bonds  # get the bonds that aren't in rings
    for bond_idx in none_ring_bonds:  # for each non-ring bond
        bgn_idx = mol.GetBondWithIdx(bond_idx).GetBeginAtomIdx()  # get the atom at each side
        end_idx = mol.GetBondWithIdx(bond_idx).GetEndAtomIdx()
        if mol.GetBondWithIdx(bond_idx).GetBondTypeAsDouble() == 1.0:  # if it is a single bond
            # get bond type as double - 1.0 for single, 1.5 for aromatic, 2.0 for double
            # if one of the atoms on the side of the bond is in the ring, add to res list
            if mol.GetAtomWithIdx(bgn_idx).IsInRing()+mol.GetAtomWithIdx(end_idx).IsInRing() == 1:
                bond = mol.GetBondWithIdx(bond_idx)
                orig_idx = bond.GetIntProp("orig_idx")
                res.append(orig_idx)
            # also if the bond atoms are in TWO DIFFERENT RINGS, add to res list
            elif not is_in_samering(bgn_idx, end_idx, bond_rings) and mol.GetAtomWithIdx(bgn_idx).IsInRing()+mol.GetAtomWithIdx(end_idx).IsInRing() == 2:
                bond = mol.GetBondWithIdx(bond_idx)
                orig_idx = bond.GetIntProp("orig_idx")
                res.append(orig_idx)
    return res



def getMaxPath(mol, path_threshold, returnMaxPath=False):
    try:
        bonds = getLinkerbond(mol)  # get the linker bonds for the molecule (join ring to linker)
        if len(bonds) == 0:
            if returnMaxPath:
                return True, None
            else:
                return True
        else:
                _fragments = Chem.FragmentOnBonds(mol, bonds)
                fragments = Chem.GetMolFrags(_fragments, asMols=True)  # split into fragments

                linkers = []
                for fragment in fragments:
                    numRings = CalcNumRings(fragment)
                    if numRings == 0:  # check if fragment is linker or ring
                        linkers.append(fragment)

                maxPath = 0
                for linker in linkers:
                    numAtoms = CalcNumHeavyAtoms(linker)
                    for i in range(numAtoms):  # look for paths up the length of number of atoms
                        # get max consecutive path length in linker
                        paths = rdmolops.FindAllPathsOfLengthN(linker, i, useBonds=True)
                        if len(paths) > 0 and i > maxPath:
                            maxPath = i

                if maxPath <= path_threshold:
                    if returnMaxPath:
                        return True, maxPath
                    else:
                        return True

                else:
                    if returnMaxPath:
                        return False, maxPath
                    else:
                        return False
    except:
        if returnMaxPath:
            return True, None
        else:
            return True


def calc_energy(mol) -> float:
    """
    Calculate energy of molecule.
    """
    mol_energy = AllChem.UFFGetMoleculeForceField(mol).CalcEnergy()
    return mol_energy


def calc_unconstrained_energy(
    og_mol, n_conf
) -> float:
    """
    Calculate average energy of multiple unconstrained conformations of a molecule.
    """
    unconstrained_energies = []
    for i in range(n_conf):  # generate conformations and calculate energy
        mol = Chem.Mol(og_mol)
        AllChem.EmbedMolecule(mol)
        AllChem.UFFOptimizeMolecule(mol)
        e = calc_energy(mol)
        unconstrained_energies.append(e)
    # calculate the average of all the energies
    avg = sum(unconstrained_energies) / len(unconstrained_energies)
    return avg


def energy_filter(mol, n_conf, energy_threshold):
    const_energy = calc_energy(mol)
    unconst_energy = calc_unconstrained_energy(mol, n_conf)

    if (const_energy / unconst_energy) >= energy_threshold:
        return False

    else:
        return True


def apply_extra_filters(names, smiles, pairs, synthons):
    mols = [Chem.MolFromSmiles(smi) for smi in smiles]

    filt_names = []
    filt_smiles = []
    filt_pairs = []
    filt_synthons = []

    for name, mol, smi, pair, synthon in zip(names, mols, smiles, pairs, synthons):
        res1 = getMaxPath(mol, 8)
        if res1:
            res2 = energy_filter(mol, 50, 7)
            if res2:
                filt_names.append(name)
                filt_smiles.append(smi)
                filt_pairs.append(pair)
                filt_synthons.append(synthon)

    return filt_names, filt_smiles, filt_pairs, filt_synthons
