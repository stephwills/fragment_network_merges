"""
Filter Fragmenstein results.
"""

import json

class FragmensteinFilter():
    """Filters the Fragmenstein results by considering how many fragments were considered
    in placing the fragment, the energy difference between the unbound and bound molecules
    and the combined RMSD between the two fragments and the placed merge."""

    def __init__(self, json_file):
        self.json_file = json_file  # the json file containing the dictionary of results
        self.data = None  # to store the dictionary
        self.result = None  # to store the filter result
    
    def get_dict(self):
        """
        Function opens the json file to load the dictionary
        """
        f = open(self.json_file)
        self.data = json.load(f)
        f.close()
    
    def filter(self):
        """
        Function filters molecules for those where both fragments were considered in its placement,
        a negative ΔΔG and combined RMSD with the fragments of < 1.5A.

        :return: returns 'pass' or 'fail'
        :rtype: string
        """
        self.get_dict()  # load the dictionary from the json
        # retrieve the energy of the bound and unbound molecules and the comRMSD
        G_bound = self.data['Energy']['ligand_ref2015']['total_score']  # energy of bound molecule
        G_unbound = self.data['Energy']['unbound_ref2015']['total_score']  # energy of unbound molecule
        comRMSD = self.data['mRMSD']  # RMSD between two fragments and merge
        regarded = 0
        # get number of molecules regarded during placement
        for rmsd in self.data['RMSDs']:
                if rmsd != None:
                    regarded += 1
        deltaG = G_bound - G_unbound  # calculate energy difference
        if regarded == 2:
            if deltaG < 0:  # keep molecules with negative ΔΔG
                if comRMSD <= 1.5:
                    self.result = 'pass'
                else:
                    self.result = 'fail'
        else:
            self.result = 'fail'
        return self.result
