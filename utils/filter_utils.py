"""
Utils functions for filtering pipeline.
"""

import os
from typing import Tuple

from pymol import cmd
from rdkit import Chem
from rdkit.Chem import AllChem, Mol, rdFMCS
from rdkit.Chem.AllChem import *


def calc_energy(mol: Mol) -> float:
    """
    Calculate energy of molecule.
    """
    mol_energy = AllChem.UFFGetMoleculeForceField(mol).CalcEnergy()
    return mol_energy


def calc_unconstrained_energy(og_mol: Mol, n_conf: int) -> float:
    """
    Calculate average energy of multiple unconstrained conformations of a molecule.
    """
    unconstrained_energies = []
    for i in range(n_conf):  # generate conformations and calculate energy
        mol = Chem.Mol(og_mol)
        AllChem.EmbedMolecule(mol, randomSeed=i*10)
        AllChem.UFFOptimizeMolecule(mol)
        e = calc_energy(mol)
        unconstrained_energies.append(e)
    # print(unconstrained_energies)
    # calculate the average of all the energies
    avg = sum(unconstrained_energies) / len(unconstrained_energies)
    # print(avg)
    return avg


def add_coordinates(fragment: Mol, substructure: Mol, atom_matches: Tuple) -> Mol:
    """
    Function to add 3D coordinates to a substructure (e.g. MCS) from the corresponding atoms from the original fragment.

    :param fragment: the original fragment (with 3D coordinates)
    :type fragment: rdkit.Chem.Mol
    :param substructure: substructure to add coordinates to
    :type substructure: rdkit.Chem.Mol
    :param atom_matches: atom numbers of the fragment that match the substructure
    :type atom_matches: tuple

    :return: substructure with coordinates added
    :rtype: rdkit.Chem.RWMol
    """
    sub_mol = Chem.RWMol(substructure)  # create editable copy of the substructure
    sub_conf = Chem.Conformer(sub_mol.GetNumAtoms())  # conformer of the substructure
    sub_matches = sub_mol.GetSubstructMatch(substructure)  # so atoms in the same order
    ref_conf = fragment.GetConformer()  # get the conformation of the actual fragment

    for i, match in enumerate(sub_matches):
        # set atom position using matching atom from fragment
        sub_conf.SetAtomPosition(match, ref_conf.GetAtomPosition(atom_matches[i]))

    sub_mol.AddConformer(sub_conf)  # add the conformation to the substructure
    return sub_mol


def remove_xe(synthon: Mol) -> Mol:
    """
    Function to remove the xenon atom from the synthon.

    :param synthon: synthon with xenon denoting attachment point
    :type synthon: rdkit.Chem.Mol

    :return: synthon with xenon removed
    :rtype: rdkit.Chem.Mol
    """
    xe = Chem.MolFromSmiles("[Xe]")
    synth = AllChem.DeleteSubstructs(synthon, xe)
    return synth


def get_mcs(full_mol: Mol, fragment: Mol) -> Mol:
    """
    Function to return the MCS between molecules. Performs partial sanitization if sanitization fails.
    """
    mcs = rdFMCS.FindMCS([full_mol, fragment], completeRingsOnly=True)
    mcs_mol = Chem.MolFromSmarts(mcs.smartsString)
    try:
        Chem.SanitizeMol(mcs_mol)
    except:
        mcs_mol.UpdatePropertyCache(strict=False)
        Chem.SanitizeMol(
            mcs_mol,
            Chem.SanitizeFlags.SANITIZE_FINDRADICALS
            | Chem.SanitizeFlags.SANITIZE_KEKULIZE
            | Chem.SanitizeFlags.SANITIZE_SETAROMATICITY
            | Chem.SanitizeFlags.SANITIZE_SETCONJUGATION
            | Chem.SanitizeFlags.SANITIZE_SETHYBRIDIZATION
            | Chem.SanitizeFlags.SANITIZE_SYMMRINGS,
            catchErrors=True,
        )
        mcs_mol.SetProp("partiallySanitized", "True")
    return mcs_mol

def sanitize(mol):
    """
    Partially sanitizes molecule if sanitization is not possible.
    """
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
        mol.SetProp("partiallySanitized", "True")


def remove_ligand(pdb_file: str) -> str:
    """
    Removes ligand from pdb file containing complex using pymol.

    :param pdb_file: pdb file of file to process
    :type pdb_file: str
    :return: name of new file
    :rtype: str
    """
    new_filename = pdb_file.replace(".pdb", "_nolig.pdb")
    if not os.path.exists(new_filename):
        cmd.reinitialize()
        cmd.load(pdb_file, "complex")
        cmd.extract("hets", "complex and HETATM")
        cmd.save(new_filename, "complex")
        # print("ligand removed")
    else:
        print("apo file exists", new_filename)
    return new_filename


def ConstrainedEmbedMatches(
    mol,
    core,
    match,
    useTethers=True,
    coreConfId=-1,
    randomseed=2342,
    getForceField=UFFGetMoleculeForceField,
    **kwargs,
):
    """
    Adapted from original RDKit function rdkit.Chem.AllChem.ConstrainedEmbed to iterate through multiple substructure
    matches
    """
    coordMap = {}
    coreConf = core.GetConformer(coreConfId)
    for i, idxI in enumerate(match):
        corePtI = coreConf.GetAtomPosition(i)
        coordMap[idxI] = corePtI

    ci = EmbedMolecule(mol, coordMap=coordMap, randomSeed=randomseed, **kwargs)
    if ci < 0:
        raise ValueError("Could not embed molecule.")

    algMap = [(j, i) for i, j in enumerate(match)]

    if not useTethers:
        # clean up the conformation
        ff = getForceField(mol, confId=0)
        for i, idxI in enumerate(match):
            for j in range(i + 1, len(match)):
                idxJ = match[j]
                d = coordMap[idxI].Distance(coordMap[idxJ])
                ff.AddDistanceConstraint(idxI, idxJ, d, d, 100.0)
        ff.Initialize()
        n = 4
        more = ff.Minimize()
        while more and n:
            more = ff.Minimize()
            n -= 1
        # rotate the embedded conformation onto the core:
        rms = AlignMol(mol, core, atomMap=algMap)
    else:
        # rotate the embedded conformation onto the core:
        rms = AlignMol(mol, core, atomMap=algMap)
        ff = getForceField(mol, confId=0)
        conf = core.GetConformer()
        for i in range(core.GetNumAtoms()):
            p = conf.GetAtomPosition(i)
            pIdx = ff.AddExtraPoint(p.x, p.y, p.z, fixed=True) - 1
            ff.AddDistanceConstraint(pIdx, match[i], 0, 0, 100.0)
        ff.Initialize()
        n = 4
        more = ff.Minimize(energyTol=1e-4, forceTol=1e-3)
        while more and n:
            more = ff.Minimize(energyTol=1e-4, forceTol=1e-3)
            n -= 1
        # realign
        rms = AlignMol(mol, core, atomMap=algMap)
    mol.SetProp("EmbedRMS", str(rms))
    return mol
