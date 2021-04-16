"""
Filters compounds using interaction fingerprints to compare interactions made with the protein with those 
made by the original fragments they were merged from.
"""

from oddt import toolkit, fingerprints

class InteractionFPFilter():

    def __init__(self, merge, fragmentA, fragmentB, protein):
        self.merge = merge
        self.fragmentA = fragmentA
        self.fragmentB = fragmentB
        self.protein = protein
    
    def get_protein(self):
        return next(toolkit.readfile('pdb', self.protein))
    
    def get_mol(self, mol):
        return next(toolkit.readfile('mol', mol))
    
    def make_fp(self, mol):
        fp = fingerprints.InteractionFingerprint(mol, self.get_protein())
        return fp
    
    def similarity_filter(self):
        protein = self.get_protein()
        merge_mol = self.get_mol(self.merge)
        fA_mol = self.get_mol(self.fragmentA)
        fB_mol = self.get_mol(self.fragmentB)
        merge_fp, fA_fp, fB_fp = self.make_fp(merge_mol), self.make_fp(fA_mol), self.make_fp(fB_mol)

        score1 = fingerprints.dice(merge_fp, fA_fp)
        score2 = fingerprints.dice(merge_fp, fB_fp)
        
        mean = (score1 + score2) / 2
        if mean >= 0.6:
            return self.merge
        else:
            return None
