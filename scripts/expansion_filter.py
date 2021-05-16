"""
Used to filter out compounds that look like elaborations rather than merges.
"""

from rdkit import Chem
from rdkit.Chem import AllChem, rdFMCS, rdShapeHelpers
from scripts.embedding_filter import add_coordinates, remove_xe

def check_for_mcs(mols):
    """Function checks if there is a MCS between fragment A, B and the merge.
    Returns the molecule if there is and returns None if not.

    :param mols: list of mols (fragments and merge)
    :type mols: list of RDKit molecules

    :return: the MCS molecule or None
    """
    mcs = rdFMCS.FindMCS(mols, completeRingsOnly=True)
    mcs_mol = Chem.MolFromSmarts(mcs.smartsString)
    smi = Chem.MolToSmiles(mcs_mol)
    if smi == '':
        return None
    else:
        mcs_mol = Chem.MolFromSmiles(Chem.MolToSmiles(mcs_mol))
        return mcs_mol

def expansion_filter(merge, fragmentA, fragmentB, synthon):
    """
    Function ensures that merges that are just expansions of one fragment are removed.
    Removes part of the merge common to both fragments (by using MCS and checking overlap).
    Then calculates the proportion of the volume of the merge that comes from expanded fragment A.
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
    synthon = remove_xe(synthon)  # remove Xe from synthon

    # in few cases there is > 1 synthon match, just pass the molecule
    if len(merge.GetSubstructMatches(synthon)) > 1:
        result = 'pass'

    else:
        # check if there is a MCS (if not returns none)
        mols = [merge, fragmentA, fragmentB]
        mcs_mol = check_for_mcs(mols)

        if mcs_mol == None:
            # if no MCS, calculate volume and ratio
            AllChem.EmbedMolecule(merge)
            merge_vol = AllChem.ComputeMolVolume(merge)
            fA_part = AllChem.DeleteSubstructs(merge, synthon)
            fA_part_vol = AllChem.ComputeMolVolume(fA_part)

            ratio = fA_part_vol / merge_vol

            if ratio >= 0.9:
                result = 'fail'
            else:
                result = 'pass'

        else:
            AllChem.EmbedMolecule(merge)  # generate conformation of the merge

            # check if there is overlap of the MCS in fragment A and B
            mcs_A = add_coordinates(fragmentA, mcs_mol)
            mcs_B = add_coordinates(fragmentB, mcs_mol)

            dist = rdShapeHelpers.ShapeProtrudeDist(mcs_A, mcs_B)

            if dist <= 0.5:  # if there is overlap, delete the substructure
                merge_no_mcs = AllChem.DeleteSubstructs(merge, mcs_mol)
                merge_vol = AllChem.ComputeMolVolume(merge_no_mcs)

                # check if synthon is still there
                if len(merge_no_mcs.GetSubstructMatches(synthon)) > 0:
                    # delete synthon and calculate the ratio
                    fA_part = AllChem.DeleteSubstructs(merge_no_mcs, synthon)
                    fA_part_vol = AllChem.ComputeMolVolume(fA_part)
                    ratio = fA_part_vol / merge_vol

                    if ratio >= 0.9:
                        result = 'fail'
                    else:
                        result = 'pass'

                else:
                    result = 'fail'

            else:   # if no overlap, calculate ratio of volumes
                merge_vol = AllChem.ComputeMolVolume(merge)
                fA_part = AllChem.DeleteSubstructs(merge, synthon)
                fA_part_vol = AllChem.ComputeMolVolume(fA_part)

                ratio = fA_part_vol / merge_vol

                if ratio >= 0.9:
                    result = 'fail'
                else:
                    result = 'pass'

    return result
