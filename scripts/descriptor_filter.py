"""
Used for filtering merges according to calculated descriptors.
"""

from rdkit import Chem
from rdkit.Chem import Crippen, rdMolDescriptors

def calculate_properties(smiles):
    """
    Function to calculate the Lipinski descriptors of the molecule.
    Counts the number of violations of Lipinski's rules.

    :param smiles: smiles of the merge
    :type smiles: smiles string

    :return: the number of Lipinski rules violated
    :rtype: int
    """
    mol = Chem.MolFromSmiles(smiles)

    # calculate the properties
    molw = rdMolDescriptors.CalcExactMolWt(mol)
    alogp = Crippen.MolLogP(mol)
    hba = rdMolDescriptors.CalcNumHBA(mol)
    hbd = rdMolDescriptors.CalcNumHBD(mol)

    # count the number of violations of Lipinski's rules
    violations = 0
    if molw > 500:
        violations += 1
    if alogp > 5:
        violations += 1
    if hba > 10:
        violations += 1
    if hbd > 5:
        violations += 1

    return violations


def descriptor_filter(smiles):
    """
    If the merge molecule violates Lipinski's rules, function returns 'fail'.

    :param smiles: smiles of the merge
    :type smiles: smiles string

    :return: returns 'pass' or 'fail'
    :rtype: string
    """
    violations = calculate_properties(smiles)
    if violations <= 1:
        result = 'pass'
    else:
        result = 'fail'
    return result
