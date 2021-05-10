"""
Used to filter out compounds that look like elaborations rather than merges.
"""

from rdkit import Chem
from rdkit.Chem import AllChem, rdFMCS, rdShapeHelpers
from scripts.embedding_filter import add_coordinates, remove_xe

def expansion_filter(merge, fragmentA, fragmentB, synthon):
    """
    Function ensures that merges that are just expansions of one fragment are removed.
    Removes part of the merge common to both fragments (by using MCS and checking overlap).
    Then calculates the proportion of the volume of the merge that comes from fragment A.
    If more than 90%, these merges are ruled out.

    :param merge: merge molecule (no 3D coordinates)
    :type merge: RDKit molecule
    :param fragmentA: fragment A molecule
    :type fragmentA: RDKit molecule
    :param fragmentB: fragment B molecule
    :type fragmentB: RDKit molecule
    :param synthon: synthon molecule (no 3D coordinates)
    :type synthon: RDKit molecule

    :return: pass or fail result
    :rtype: string
    """
    synthon = remove_xe(synthon)  # remove xenon atom from synthon

    # find common substructure between the fragments
    mcs = rdFMCS.FindMCS([synthon, fragmentA, fragmentB])
    mcs_mol = Chem.MolFromSmarts(mcs.smartsString)

    if len(fragmentA.GetSubstructMatches(mcs_mol)) <= 1:
        # check the the common substructure is actually overlapping
        # get the substructure (with coordinates) from fragment A and B by removing everything else
        substructA = Chem.ReplaceSidechains(fragmentA, mcs_mol)
        substructB = Chem.ReplaceSidechains(fragmentB, mcs_mol)
        # check they are overlapping with protrusion distance
        dist = rdShapeHelpers.ShapeProtrudeDist(substructA, substructB)

        if dist <= 0.5:  # if overlapping by at least 50%
            # remove the common part of the merge
            merge_no_mcs = AllChem.DeleteSubstructs(merge, mcs_mol)
            # create random conformer (for the computation of mol volume)
            AllChem.EmbedMolecule(merge_no_mcs)
            merge_no_mcs_vol = AllChem.ComputeMolVolume(merge_no_mcs)

            # isolate the part of the molecule coming from fragment A
            # (after removing parts common to both molecules)
            fA_part = AllChem.DeleteSubstructs(merge, synthon)
            fA_part = add_coordinates(merge_no_mcs, fA_part)
            fA_part_vol = AllChem.ComputeMolVolume(fA_part)  # calculate volume

            # calculate how much fragment A contributes to the total volume
            ratio = fA_part_vol / merge_no_mcs_vol
            if ratio >= 0.9:  # if more than 90%, do not use
                result = 'fail'
            else:
                result = 'pass'
        else:
            result = 'pass'
    
    else:
        result = 'pass'
    
    return result

def simple_filter(merge, fragmentA, fragmentB, synthon):
    """
    Filter to calculate proportion of molecule that came from fragment A.
    """
    synthon = remove_xe(synthon)
    merge_mol = Chem.Mol(merge)
    AllChem.EmbedMolecule(merge_mol)
    
    merge_vol = AllChem.ComputeMolVolume(merge_mol)
    fA_part = AllChem.DeleteSubstructs(merge_mol, synthon)
    fA_part_vol = AllChem.ComputeMolVolume(fA_part)

    ratio = fA_part_vol / merge_vol
    if ratio >= 0.9:
        result = 'fail'
    else:
        result = 'pass'
    
    return result
