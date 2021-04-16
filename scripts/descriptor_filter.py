"""
Used for filtering merges according to calculated descriptors.
"""

from rdkit import Chem
from rdkit.Chem import Crippen, rdMolDescriptors

class DescriptorFilter():
    """Filters merges according to whether calculated descriptors meet some defined criteria."""

    def __init__(self, smiles):
        self.smiles = smiles
        self.result = None

    def calculate_properties(self):
        """
        Function to calculate various QED properties of the compounds.
        Checks the properties against some defined criteria and counts the number of 'violations'.

        :return: the number of 'violations' of the descriptor criteria
        :rtype: int
        """
        mol = Chem.MolFromSmiles(self.smiles)
        # calculate the properties
        molw = rdMolDescriptors.CalcExactMolWt(mol)
        alogp = Crippen.MolLogP(mol)
        hba = rdMolDescriptors.CalcNumHBA(mol)
        hbd = rdMolDescriptors.CalcNumHBD(mol)
        psa = rdMolDescriptors.CalcTPSA(mol)
        rotat = rdMolDescriptors.CalcNumRotatableBonds(mol)
        rings = rdMolDescriptors.CalcNumAromaticRings(mol)
        violations = 0
        # count the number of violations of some defined criteria (this can be modified)
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

    def filter(self):
        """
        Function filters the smiles according to the number of descriptor violations.
        If more than two violations, the function returns None.

        :return: returns the smiles if the molecule passes; otherwise None
        :rtype: smiles string or None
        """
        violations = self.calculate_properties()
        if violations <= 2:
            self.result = 'pass'
        else:
            self.result = 'fail'
        return self.result
