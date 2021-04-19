"""
Used for filtering merges according to calculated descriptors.
"""

from rdkit import Chem
from rdkit.Chem import Crippen, rdMolDescriptors

def calculate_properties(smiles):
    """
    Function to calculate various QED properties of the compounds.
    Checks the properties against some defined criteria and counts the number of 'violations'.

    :param smiles: smiles of the merge
    :type smiles: smiles string

    :return: the number of 'violations' of the descriptor criteria
    :rtype: int
    """
    mol = Chem.MolFromSmiles(smiles)
    # calculate the properties
    molw = rdMolDescriptors.CalcExactMolWt(mol)
    alogp = Crippen.MolLogP(mol)
    hba = rdMolDescriptors.CalcNumHBA(mol)
    hbd = rdMolDescriptors.CalcNumHBD(mol)
    psa = rdMolDescriptors.CalcTPSA(mol)
    rotat = rdMolDescriptors.CalcNumRotatableBonds(mol)
    rings = rdMolDescriptors.CalcNumAromaticRings(mol)
    violations = 0
    # count the number of violations of the defined criteria (this can be modified)
    if molw > 350:
        violations += 1
    if alogp > 3:
        violations += 1
    if hba > 3:
        violations += 1
    if hbd > 3:
        violations += 1
    if psa > 100:
        violations += 1
    if rotat > 8:
        violations += 1
    if rings > 2:
        violations += 1
    return violations

def descriptor_filter(smiles):
    """
    Function filters the smiles according to the number of descriptor violations.
    If more than two violations, the function returns 'fail'.

    :param smiles: smiles of the merge
    :type smiles: smiles string

    :return: returns 'pass' or 'fail'
    :rtype: string
    """
    violations = calculate_properties(smiles)
    if violations <= 2:
        result = 'pass'
    else:
        result = 'fail'
    return result
