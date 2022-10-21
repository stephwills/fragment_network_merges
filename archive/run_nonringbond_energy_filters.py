"""Extra filters"""

import os
from rdkit import Chem, RDLogger
from rdkit.Chem.Scaffolds import MurckoScaffold
from rdkit.Chem.rdMolDescriptors import CalcNumRings
from rdkit.Chem import rdmolops, AllChem
from filter.config_filter import config_filter
from joblib import Parallel, delayed
from tqdm import tqdm

# code from https://iwatobipen.wordpress.com/2020/01/23/cut-molecule-to-ring-and-linker-with-rdkit-rdkit-chemoinformatics-memo/
RDLogger.DisableLog("rdApp.*")


def is_in_samering(idx1, idx2, bond_rings):
    for bond_ring in bond_rings:
        if idx1 in bond_ring and idx2 in bond_ring:
            return True
    return False


def getLinkerbond(mol, useScaffold=True):
    res = []
    for atom in mol.GetAtoms():
        atom.SetIntProp("orig_idx", atom.GetIdx())  # set original atom index as property
    for bond in mol.GetBonds():
        bond.SetIntProp("orig_idx", bond.GetIdx())  # set original bond index as property

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


def get_fragments(mol, bonds):
    _fragments = Chem.FragmentOnBonds(mol, bonds)
    fragments = Chem.GetMolFrags(_fragments, asMols=True)
    linkers = []
    sidechains = []

    for fragment in fragments:
        # print(fragments)
        numRings = CalcNumRings(fragment)
        if numRings == 0:
            smiles = Chem.MolToSmiles(fragment)
            # print(smiles)
            if smiles.count('*') == 1:
                sidechains.append(fragment)
            else:
                linkers.append(fragment)

    return linkers, sidechains

def max_path(fragment):
    maxPath = 0
    numAtoms = fragment.GetNumHeavyAtoms()
    for i in range(numAtoms):  # look for paths up the length of number of atoms
        # get max consecutive path length in linker
        paths = rdmolops.FindAllPathsOfLengthN(fragment, i, useBonds=True)
        if len(paths) > 0:
            maxPath = i

    return maxPath


def getMaxPath(mol,
               linkerPathThreshold=config_filter.LINKER_PATH_THRESHOLD,
               sidechainPathThreshold=config_filter.SIDECHAIN_PATH_THRESHOLD,
               returnMaxPath=False,
               useScaffold=False):
    try:
        bonds = getLinkerbond(mol, useScaffold)  # get the linker bonds for the molecule (join ring to linker)
        if len(bonds) == 0:
            if returnMaxPath:
                return True, None
            else:
                return True
        else:
            linkers, sidechains = get_fragments(mol, bonds)

            maxLinkerPath = 0
            maxSidechainPath = 0

            if len(linkers) > 0:
                for linker in linkers:
                    numLinkerAtoms = linker.GetNumHeavyAtoms()
                    for i in range(numLinkerAtoms):  # look for paths up the length of number of atoms
                        # get max consecutive path length in linker
                        linkerPaths = rdmolops.FindAllPathsOfLengthN(linker, i, useBonds=True)
                        if len(linkerPaths) > 0 and i > maxLinkerPath:
                            maxLinkerPath = i

            if len(sidechains) > 0:
                for sidechain in sidechains:
                    numSidechainAtoms = sidechain.GetNumHeavyAtoms()
                    for i in range(numSidechainAtoms):  # look for paths up the length of number of atoms
                        # get max consecutive path length in sidechain
                        sidechainPaths = rdmolops.FindAllPathsOfLengthN(sidechain, i, useBonds=True)
                        if len(sidechainPaths) > 0 and i > maxSidechainPath:
                            maxSidechainPath = i

            if maxLinkerPath < linkerPathThreshold and maxSidechainPath < sidechainPathThreshold:
                if returnMaxPath:
                    return True, maxLinkerPath, maxSidechainPath
                else:
                    return True

            else:
                if returnMaxPath:
                    return False, maxLinkerPath, maxSidechainPath
                else:
                    return False
    except:
        if returnMaxPath:
            return True, None, None
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


def get_energy_ratio(mol, n_conf, energy_threshold=config_filter.ENERGY_THRESHOLD):
    const_energy = calc_energy(mol)
    unconst_energy = calc_unconstrained_energy(mol, n_conf)
    return const_energy, unconst_energy, (const_energy / unconst_energy)


def energy_filter(mol, n_conf, energy_threshold=config_filter.ENERGY_THRESHOLD):
    const_energy = calc_energy(mol)
    unconst_energy = calc_unconstrained_energy(mol, n_conf)

    if (const_energy / unconst_energy) >= energy_threshold:
        return False

    else:
        return True


def extra_filters(mol):
    res1 = getMaxPath(mol, useScaffold=False)
    if res1:
        try:
            res2 = energy_filter(mol, 50, 7)
        except:
            return False
        if res2:
            return True
        else:
            return False
    else:
        return False


def apply_maxpath_filter(smiles):
    mols = [Chem.MolFromSmiles(s) for s in smiles]
    results = Parallel(n_jobs=os.cpu_count(), backend='multiprocessing')(
        delayed(getMaxPath)(mol, useScaffold=False) for mol in tqdm(mols))
    return results


def apply_extra_filters(names, smiles, pairs, synthons, ifp_scores, sucos_scores, mol_files, sim_search=False):
    n_original = len(smiles)
    mols = [Chem.rdmolfiles.MolFromMolFile(i) for i in mol_files]

    filt_names = []
    filt_smiles = []
    filt_pairs = []
    if not sim_search:
        filt_synthons = []
    filt_ifp_scores = []
    filt_sucos_scores = []

    results = Parallel(n_jobs=os.cpu_count(), backend='multiprocessing')(delayed(extra_filters)(mol) for mol in tqdm(mols))
    if not sim_search:
        for name, smi, pair, synthon, ifp_score, sucos_score, res in zip(names, smiles, pairs, synthons, ifp_scores, sucos_scores, results):
            if res:
                filt_names.append(name)
                filt_smiles.append(smi)
                filt_pairs.append(pair)
                filt_synthons.append(synthon)
                filt_ifp_scores.append(ifp_score)
                filt_sucos_scores.append(sucos_score)
        n_filtered = len(filt_smiles)
        print(f"{n_filtered}/{n_original} smiles remaining")
        return filt_names, filt_smiles, filt_pairs, filt_synthons, filt_ifp_scores, filt_sucos_scores

    else:
        for name, smi, pair, ifp_score, sucos_score, res in zip(names, smiles, pairs, ifp_scores, sucos_scores, results):
            if res:
                filt_names.append(name)
                filt_smiles.append(smi)
                filt_pairs.append(pair)
                filt_ifp_scores.append(ifp_score)
                filt_sucos_scores.append(sucos_score)

        n_filtered = len(filt_smiles)

        print(f"{n_filtered}/{n_original} smiles remaining")
        return filt_names, filt_smiles, filt_pairs, filt_ifp_scores, filt_sucos_scores
