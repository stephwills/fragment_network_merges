"""
Used to generate merges between two fragments using the fragment network.
Uses code from https://github.com/tdudgeon/fragment-network-merges.
"""

import itertools
import getpass
from neo4j import GraphDatabase
from rdkit import Chem

try:
    password
except NameError:
    password = getpass.getpass()

driver = GraphDatabase.driver("bolt://localhost:7687", auth=("swills", password))

def get_combinations(fragments):
    """
    Enumerate all possible combinations of fragments for merging.

    :return: list of possible combinations of fragments
    :rtype: list of tuples
    """
    permutations = list(itertools.permutations(fragments, 2))
    return permutations

def mols_from_smiles(list_smiles):
    """
    Generates list of RDKit molecules from list of smiles.

    :param list_smiles: list of smiles strings
    :type list_smiles: list

    :return: list of rdkit molecules
    :rtype: list
    """
    mols = [Chem.MolFromSmiles(smiles) for smiles in list_smiles]
    return mols

def find_molecule_node(tx, smiles):
    """
    Finds node in the fragment network.

    :param tx: transaction in which query is run
    :type tx: neo4j transaction
    :param smiles: smiles of the fragment
    :type smiles: string

    :return: molecule node
    :rtype: node
    """
    for record in tx.run('MATCH (m:F2 {smiles: $smiles}) RETURN m', smiles=smiles):
        node = record['m']
        return node

def check_for_nodes(fragments):
    """
    Checks if nodes are present.

    :param fragments: list of fragment smiles
    :type fragments: list
    """
    number_nodes = 0
    number_fragments = len(fragments)
    with driver.session() as session:
        for fragment in fragments:
            mol = session.read_transaction(find_molecule_node, fragment)  # run neo4j query
            if mol:
                number_nodes += 1
        print(f'{number_nodes} of {number_fragments} fragments present in network')

def filter_for_nodes(fragments):
    """
    Checks if nodes are present and filters fragments for those that are in network.

    :param fragments: list of fragment smiles
    :type fragments: list

    :return: list of available fragments
    :rtype: list of smiles strings
    """
    filtered_fragments = []
    with driver.session() as session:
        for fragment in fragments:
            mol = session.read_transaction(find_molecule_node, fragment)  # run neo4j query
            if mol:
                filtered_fragments.append(fragment)
    return filtered_fragments

def add_required_synthons(labels, synthon):
    """
    Checks that the synthons have a single attachment point and single component.

    :param labels: set containing the synthons
    :type labels: set
    :param synthon: smiles string containing synthon with Xe atom
    :type synthon: string
    """
    if synthon.count('[Xe]') > 0:
        labels.add(synthon)

def find_synthons(tx, smiles):
    """Query for all child fragments (recursive).
    Extract the label property of each edge and collect a set of SMILES that match our needs.
    """
    labels = set()
    for record in tx.run('MATCH (fa:F2 {smiles: $smiles})-[e:FRAG*]->(f:F2) RETURN e',
                         smiles=smiles):
        edges = record['e']
        for edge in edges:
            s = edge['label']
            tokens = s.split('|')
            #print('Found', tokens[1], tokens[4])
            add_required_synthons(labels, tokens[1])
            add_required_synthons(labels, tokens[4])
    return list(labels)

def get_synthons(smiles):
    """
    Extract the synthons from the database.

    :param smiles: smiles of the fragment
    :type smiles: string

    :return: synthons
    :rtype: string
    """
    with driver.session() as session:
        synthons = session.read_transaction(find_synthons, smiles)
        print(f'Found {len(synthons)} synthons')
        return synthons

def find_expansions(tx, smiles, synthon):
    """
    Expand fragment 'A' using the synthons generated from fragment 'B'
    (or vice versa).

    :param smiles: smiles of the fragment to expand
    :type smiles: string
    :param synthon: synthon of the generated synthon
    :type synthon: string

    :return: expansions
    :rtype: set
    """
    query = ("MATCH (fa:F2 {smiles: $smiles})"
             "-[:FRAG*0..2]-(:F2)"
             "<-[e:FRAG]-(c:Mol) WHERE"
             " c.hac > 15 AND"
             " (split(e.label, '|')[1] = $synthon OR split(e.label, '|')[4] = $synthon)"
             " RETURN DISTINCT c")
    expansions = set()
    for record in tx.run(query, smiles=smiles, synthon=synthon):
        node = record['c']
        expansions.add(node['smiles'])
    return expansions

def get_expansions(fragment, synthons):
    """
    Get the expansions of when fragment using a set of synthons from the other fragment.
    Returns a dictionary where the keys are the synthon used in the expansion, and the values
    are a list of the expansions.

    :param fragment: smiles string of the fragment to expand
    :type fragment: string
    :param synthons: list of synthons smiles
    :type synthons: list

    :return: dictionary containing expansions
    :rtype: dictionary
    """
    all_expansions = {}
    with driver.session() as session:
        number = 0
        total_expansions = 0
        expanded_synthons = 0
        for synthon in synthons:
            expansions = session.read_transaction(find_expansions, fragment, synthon)
            all_expansions[synthon] = expansions
            print(f'Synthon {number}: found {len(expansions)} expansions')
            number += 1
            total_expansions += len(expansions)
            if expansions:  # record if the synthon led to expansions
                expanded_synthons += 1
    print(f'{total_expansions} expansions from {expanded_synthons} out of {len(synthons)} synthons')
    print('')
    return all_expansions
