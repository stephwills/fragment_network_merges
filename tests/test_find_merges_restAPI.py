
from rdkit import Chem

smi = Chem.MolToSmiles(Chem.MolFromSmiles("CC=1C=CC(CS(=O)(=O)N)=CC1F"))
print(smi)
merger = MergerFinder_restAPI()
results = merger.check_for_nodes([smi])
print(results)

results = merger.get_synthons(smi)
print(results)

smi2 = Chem.MolToSmiles(Chem.MolFromSmiles("CS(=O)(=O)NCCC=1C=CC=CC1"))
results = merger.get_expansions([smi, smi2], ["x0034_0B", "x0176_0B"], "nsp13", None)
print(results)