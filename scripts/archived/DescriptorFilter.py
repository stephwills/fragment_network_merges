"""
Used for filtering fragment merges.
"""

from rdkit import Chem
from rdkit.Chem import AllChem, Crippen, rdFMCS, rdForceFieldHelpers, rdMolDescriptors, PandasTools
import numpy as np

# add documentation

class DescriptorFilter():

    def __init__(self, merges):
        self.merges = merges

    def calculate_properties(self, smiles):
        mol = Chem.MolFromSmiles(smiles)
        mw = rdMolDescriptors.CalcExactMolWt(mol)
        alogp = Crippen.MolLogP(mol)
        hba = rdMolDescriptors.CalcNumHBA(mol)
        hbd = rdMolDescriptors.CalcNumHBD(mol)
        psa = rdMolDescriptors.CalcTPSA(mol)
        rotat = rdMolDescriptors.CalcNumRotatableBonds(mol)
        rings = rdMolDescriptors.CalcNumAromaticRings(mol)
        violations = 0
        if mw > 350: violations += 1
        if alogp > 3: violations += 1
        if hba > 3: violations += 1
        if hbd > 3: violations += 1
        if psa > 100: violations += 1
        if rotat > 8: violations += 1
        if rings > 2: violations += 1
        return violations
 
    def filter(self):
        for smi in self.merges:
            violations = self.calculate_properties(smi)
            if violations > 2:
                self.merges.remove(smi)
        
        return self.merges
  
    def filter_bool(self):
        boolean = []
        for smi in self.merges:
            violations = self.calculate_properties(smi)
            if violations > 2:
                boolean.append(True)
            else:
                boolean.append(False)
        return boolean

class DescriptorFilterDf():

    def __init__(self, df, smiCol='Smiles', synthCol='Synthon', keepProperties=False):
        self.df = df
        self.smiCol = smiCol
        self.synthCol = synthCol
        self.keepProperties = keepProperties
    
    def add_mols(self):
        PandasTools.AddMoleculeColumnToFrame(self.df, smilesCol=self.smiCol, molCol='Molecules')
        PandasTools.AddMoleculeColumnToFrame(self.df, smilesCol=self.synthCol, molCol='Synthon molecules')

    def calc_properties(self):
        self.add_mols()
        self.df['Molecular weight'] = self.df['Molecules'].apply(rdMolDescriptors.CalcExactMolWt)
        self.df['AlogP'] = self.df['Molecules'].apply(Crippen.MolLogP)
        self.df['HBA'] = self.df['Molecules'].apply(rdMolDescriptors.CalcNumHBA)
        self.df['HBD'] = self.df['Molecules'].apply(rdMolDescriptors.CalcNumHBD)
        self.df['PSA'] =  self.df['Molecules'].apply(rdMolDescriptors.CalcTPSA)
        self.df['Rotatable bonds'] = self.df['Molecules'].apply(rdMolDescriptors.CalcNumRotatableBonds)
        self.df['Aromatic rings'] = self.df['Molecules'].apply(rdMolDescriptors.CalcNumAromaticRings)
    
    def calc_violations(self):
        self.calc_properties()
        num_violations = []
        for _, row in self.df.iterrows():
            violations = 0
            if row['Molecular weight'] > 350: violations += 1
            if row['AlogP'] > 3: violations += 1
            if row['HBA'] > 3: violations += 1
            if row['HBD'] > 3: violations += 1
            if row['PSA'] > 100: violations += 1
            if row['Rotatable bonds'] > 8: violations += 1
            if row['Aromatic rings'] > 2: violations += 1
            num_violations.append(violations)
        self.df['Violations'] = num_violations
    
    def filter_df(self):
        self.calc_violations()
        no_compounds = len(self.df)
        self.df = self.df[self.df.Violations <= 2].reset_index(drop=True)
        removed = no_compounds - len(self.df)
        print(f'Removed {removed} compounds. {len(self.df)} remaining.')
        if self.keepProperties == False:
            propertyCols = ['Violations', 'Molecular weight', 'AlogP', 'HBA', 'HBD', 'PSA', 'Rotatable bonds', 'Aromatic rings']
            self.df = self.df.drop(propertyCols, axis=1)

        return self.df
