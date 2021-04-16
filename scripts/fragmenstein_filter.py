"""
Filter Fragmenstein results.
"""

import json

class FragmensteinFilter():

    def __init__(self, json_file):
        self.json_file = json_file
        self.data = None
        self.result = None
    
    def get_dict(self):
        f = open(self.json_file)
        self.data = json.load(f)
        f.close()
    
    def filter(self):
        self.get_dict()
        G_bound = self.data['Energy']['ligand_ref2015']['total_score']
        G_unbound = self.data['Energy']['unbound_ref2015']['total_score']
        comRMSD = self.data['mRMSD']
        regarded = 0
        for rmsd in self.data['RMSDs']:
                if rmsd != None:
                    regarded += 1
        deltaG = G_bound - G_unbound
        if regarded == 2:
            if deltaG < 0:
                if comRMSD <= 1.5:
                    self.result = 'pass'
                else:
                    self.result = 'fail'
        else:
            self.result = 'fail'
        
        return self.result
