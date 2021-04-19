"""
Filter Fragmenstein results.
"""

import json
    
def get_dict(json_file):
    """
    Function opens the json file to load the dictionary

    :param json_file: json containing the dictionary
    :type json_file: .json
    """
    f = open(json_file)
    data = json.load(f)
    f.close()
    return data

def fragmenstein_filter(json_file):
    """
    Function filters molecules for those where both fragments were considered in its placement,
    a negative ΔΔG and combined RMSD with the fragments of < 1.5A.

    :param json_file: json containing the dictionary
    :type json_file: .json

    :return: returns 'pass' or 'fail'
    :rtype: string
    """
    data = get_dict(json_file) # load the dictionary from the json
    # retrieve the energy of the bound and unbound molecules and the comRMSD
    G_bound = data['Energy']['ligand_ref2015']['total_score']  # energy of bound molecule
    G_unbound = data['Energy']['unbound_ref2015']['total_score']  # energy of unbound molecule
    comRMSD = data['mRMSD']  # RMSD between two fragments and merge
    regarded = 0
    # get number of molecules regarded during placement
    for rmsd in data['RMSDs']:
            if rmsd != None:
                regarded += 1
    deltaG = G_bound - G_unbound  # calculate energy difference
    if regarded == 2:
        if deltaG < 0:  # keep molecules with negative ΔΔG
            if comRMSD <= 1.5:
                result = 'pass'
            else:
                result = 'fail'
    else:
        result = 'fail'
    return result
