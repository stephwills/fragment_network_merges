"""
Used to generate merges between two fragments using the fragment network.
Uses code from https://github.com/tdudgeon/fragment-network-merges.
"""

import os
import itertools
import getpass
import json
from neo4j import GraphDatabase
from rdkit import Chem
from rdkit.Chem import rdShapeHelpers
from scripts.embedding_filter import add_coordinates, remove_xe
from scripts.preprocessing import get_mol

password = getpass.getpass()
driver = GraphDatabase.driver("bolt://localhost:7687", auth=("swills", password))

# functions for checking the nodes and filtering for fragments that exist as nodes
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
    Prints statement with how many of the fragments are present in the network.

    :param fragments: list of smiles of the fragments
    :type fragments: list of strings
    """
    number_nodes = 0
    number_fragments = len(fragments)
    with driver.session() as session:
        for fragment in fragments:
            mol = session.read_transaction(find_molecule_node, fragment)  # run neo4j query
            if mol:
                number_nodes += 1
        print(f'{number_nodes} of {number_fragments} fragments present in network')

def filter_for_nodes(fragments, names):
    """
    Checks if nodes are present and filters fragments for those that are in network.

    :param fragments: list of fragment smiles
    :type fragments: list of strings
    :param names: list of fragment names
    :type names: list of strings

    :return: list of available fragments
    :rtype: list of smiles strings
    :return: list of available fragments (names)
    :rtype: list of strings
    """
    removed = 0
    with driver.session() as session:
        for i, fragment in enumerate(fragments):
            mol = session.read_transaction(find_molecule_node, fragment)  # run neo4j query
            if not mol:
                fragments.pop(i)
                names.pop(i)
                removed += 1
    print(f'{removed} fragments removed from list. {len(fragments)} fragments remaining.')
    return fragments, names

def get_combinations(fragments, names):
    """
    Enumerate all possible combinations of fragments for merging.

    :param fragments: list of fragment smiles
    :type fragments: list of strings
    :param names: list of fragment names
    :type names: list of strings

    :return: list of possible combinations of fragments
    :rtype: list of tuples
    :return: list of possible combinations of fragments (names)
    :rtype: list of tuples
    """
    fragment_pairs = list(itertools.permutations(fragments, 2))
    name_pairs = list(itertools.permutations(names, 2))
    return fragment_pairs, name_pairs

# code for generating the list of fragment merges
def add_required_synthons(labels, synthon):
    """
    Checks that the synthons have a single attachment point and single component.

    :param labels: set containing the synthons
    :type labels: set
    :param synthon: smiles string containing synthon with Xe atom
    :type synthon: string
    """
    if synthon.count('[Xe]') == 1:
        labels.add(synthon)

def find_synthons(tx, smiles):
    """
    Query for all child fragments (recursive).
    Extract the label property of each edge and collect a set of SMILES that match our needs.
    """
    labels = set()
    for record in tx.run('MATCH (fa:F2 {smiles: $smiles})-[e:FRAG*]->(f:F2) RETURN e',
                            smiles=smiles):
        edges = record['e']
        for edge in edges:
            s = edge['label']
            tokens = s.split('|')
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
    Expand fragment 'A' using the synthons generated from fragment 'B' using a neo4j
    query. Query limited to compounds available from vendors, with HAC > 15
    and a maximum of 2 hops away.

    :param smiles: smiles of the fragment to expand
    :type smiles: string
    :param synthon: synthon of the generated synthon
    :type synthon: string

    :return: expansions
    :rtype: set
    """
    # query = ("MATCH (fa:F2 {smiles: $smiles})"
    #             "-[:FRAG*0..2]-(:F2)"
    #             "<-[e:FRAG]-(c:Mol) WHERE"
    #             " c.hac > 15 AND"
    #             " (split(e.label, '|')[1] = $synthon OR split(e.label, '|')[4] = $synthon)"
    #             " RETURN DISTINCT c")
    query = ("MATCH (fa:F2 {smiles: $smiles})"
                "-[:FRAG]->(:F2)"
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

def filter_synthons(synthon):
    """
    Counts the number of carbons in the synthon. Keep synthons that have at least three carbons.

    :param synthon: synthon of interest
    :type synthon: smiles string

    :return: synthon smiles or None
    :rtype: smiles string or None
    """
    carbon = Chem.MolFromSmarts('[#6]')
    synthon_mol = Chem.MolFromSmiles(synthon)
    num_carbons = len(synthon_mol.GetSubstructMatches(carbon))
    if num_carbons >= 3:
        return synthon

def substructure_check(synthon, fragmentA, fragmentB):
    """
    Checks if the synthon is already present in fragment A and in an overlapping position.
    These synthons result in elaborations of fragment A rather than a merge.

    :param synthon: synthon of interest
    :type synthon: smiles string
    :param fragment A: fragment undergoing expansions
    :type fragment A: RDKit molecule

    :return: synthon smiles or None
    :rtype: smiles string or None
    """
    # add coordinates from the fragments to the synthon structure
    synthon_A = Chem.MolFromSmiles(synthon)
    synthon_A = remove_xe(synthon_A)
    synthon_B = Chem.Mol(synthon_A)

    try:
        fA_match = add_coordinates(fragmentA, synthon_A)
        fB_match = add_coordinates(fragmentB, synthon_B)

        distance = rdShapeHelpers.ShapeProtrudeDist(fA_match, fB_match)
        # if they overlap by more than 50%, then remove
        if distance >= 0.8:
            return synthon
        else:
            return None

    except:
        return synthon

def get_expansions(fragments, names, target, output_dir):
    """
    Function executes the whole process, generating synthons for fragment B and using them to
    generate expansions of fragment A. Returns a dictionary containing all the synthons as keys,
    and the smiles of the merges they were used to generate.
    The dictionary is then saved as a json file.

    :return: dictionary of merges
    :rtype: dictionary
    """
    # get fragment A and B from the tuple
    fragmentA, fragmentB = fragments[0], fragments[1]
    nameA, nameB = names[0], names[1]
    molA, molB = get_mol(target, nameA), get_mol(target, nameB)
    print(f'Expanding fragment A: {nameA} with synthons of fragment B: {nameB}')

    # generate the synthons from fragment B
    unfiltered_synthons = get_synthons(fragmentB)

    # filter synthons for those with <3 carbons
    synthons_3c = [filter_synthons(syn) for syn in unfiltered_synthons]
    # remove None values from list
    synthons_3c = list(filter(None, synthons_3c))

    # filter synthons for those already in fragment A
    synthons = [substructure_check(syn, molA, molB) for syn in synthons_3c]
    # remove None values from list
    synthons = list(filter(None, synthons))
    print(f'{len(synthons)} synthons remaining after filtering')

    # create empty dictionary to store results
    all_expansions = {}
    with driver.session() as session:
        number = 0
        total_expansions = 0
        expanded_synthons = 0
        for synthon in synthons:  # expand fragment A using each synthon
            print(f'Running synthon {number}')
            expansions = session.read_transaction(find_expansions, fragmentA, synthon)
            all_expansions[synthon] = list(expansions)  # store in dictionary with the synthon as key
            print(f'Synthon {number}: found {len(expansions)} expansions')
            number += 1
            total_expansions += len(expansions)
            if expansions:  # record if the synthon led to expansions
                expanded_synthons += 1
    print(f'{total_expansions} expansions from {expanded_synthons} out of {len(synthons)} synthons')

    # save as json file
    filename = nameA + '_' + nameB + '.json'
    filepath = os.path.join(output_dir, filename)
    with open(filepath, 'w') as f:
        json.dump(all_expansions, f)
