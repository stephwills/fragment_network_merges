"""
Filters compounds using interaction fingerprints to compare interactions made with the protein with those 
made by the original fragments they were merged from.
"""

from oddt import toolkit, fingerprints

# check fp calculation; change to Tversky similarity?
    
def get_protein(protein):
    """
    Function loads the protein from the pdb file.

    :param protein: the protein file
    :type protein: pdb file

    :return: protein
    :rtype: ODDT protein
    """
    return next(toolkit.readfile('pdb', protein))

def get_mol(mol):
    """
    Function loads the molecule from the mol file.

    :param mol: mol to read (fragment or placed merge)
    :type mol: mol file

    :return: molecule
    :rtype: ODDT molecule
    """
    return next(toolkit.readfile('mol', mol))

def make_fp(protein, mol):
    """
    Function creates interaction fingerprint between the molecule and protein.

    :param protein: protein loaded from file
    :type protein: ODDT protein
    :param mol: mol loaded from file
    :type mol: ODDT mol

    :return: fingerprint
    :rtype: list
    """
    fp = fingerprints.InteractionFingerprint(mol, protein)
    return fp

def similarity_filter(merge, fragmentA, fragmentB, prot):
    """
    Function compares the interaction fingerprints between the fragment and the protein,
    and the merge and the protein. Calculates similarity (using dice). Keeps molecules
    with similarity coefficient of greater than 0.6.

    :param merge: the file of the merge to filter
    :type merge: mol file
    :param fragmentA: the fragment A file
    :type fragmentA: mol file
    :param fragmentB: the fragment B file
    :type fragmentB: mol file
    :param prot: protein file
    :type prot: pdb file

    :return: returns 'pass' or 'fail'
    :rtype: string
    """
    # load the molecules
    protein = get_protein(prot)
    merge_mol = get_mol(merge)
    fA_mol = get_mol(fragmentA)
    fB_mol = get_mol(fragmentB)
    # create all the fingerprints
    merge_fp, fA_fp, fB_fp = make_fp(protein, merge_mol), make_fp(protein, fA_mol), make_fp(protein, fB_mol)
    # calculate the dice coefficient between the fp for the merge and the two fragments
    score1 = fingerprints.dice(merge_fp, fA_fp)
    score2 = fingerprints.dice(merge_fp, fB_fp)
    # calculate the average between the two scores
    mean = (score1 + score2) / 2
    if mean >= 0.6:  # keep molecules with similarity > 0.6
        result = 'pass'
    else:
        result = 'fail'
    return result
