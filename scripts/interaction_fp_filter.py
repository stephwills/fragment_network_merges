"""
Filters compounds using interaction fingerprints to compare interactions made with the protein with those 
made by the original fragments they were merged from.
"""

from oddt import toolkit, fingerprints

# check fp calculation; change to Tversky similarity?
class InteractionFPFilter():
    """Filters molecules for those that maintain the interactions made by the fragments
    with the residues"""

    def __init__(self, merge, fragmentA, fragmentB, protein):
        self.merge = merge
        self.fragmentA = fragmentA
        self.fragmentB = fragmentB
        self.protein = protein
        self.result = None
    
    def get_protein(self):
        """
        Function loads the protein from the pdb file.

        :return: protein
        :rtype: ODDT protein
        """
        return next(toolkit.readfile('pdb', self.protein))
    
    def get_mol(self, mol):
        """
        Function loads the molecule from the mol file.

        :return: molecule
        :rtype: ODDT molecule
        """
        return next(toolkit.readfile('mol', mol))
    
    def make_fp(self, mol):
        """
        Function creates interaction fingerprint between the molecule and protein.

        :return: fingerprint
        :rtype: list
        """
        fp = fingerprints.InteractionFingerprint(mol, self.get_protein())
        return fp
    
    def similarity_filter(self):
        """
        Function compares the interaction fingerprints between the fragment and the protein,
        and the merge and the protein. Calculates similarity (using dice). Keeps molecules
        with similarity coefficient of greater than 0.6.

        :return: returns 'pass' or 'fail'
        :rtype: string
        """
        # load the molecules
        protein = self.get_protein()
        merge_mol = self.get_mol(self.merge)
        fA_mol = self.get_mol(self.fragmentA)
        fB_mol = self.get_mol(self.fragmentB)
        # create all the fingerprints
        merge_fp, fA_fp, fB_fp = self.make_fp(merge_mol), self.make_fp(fA_mol), self.make_fp(fB_mol)
        # calculate the dice coefficient between the fp for the merge and the two fragments
        score1 = fingerprints.dice(merge_fp, fA_fp)
        score2 = fingerprints.dice(merge_fp, fB_fp)
        # calculate the average between the two scores
        mean = (score1 + score2) / 2
        if mean >= 0.6:  # keep molecules with similarity > 0.6
            self.result = 'pass'
        else:
            self.result = 'fail'
        return self.result
