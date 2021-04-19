"""
Used to generate merges between two fragments using the fragment network.
Uses code from https://github.com/tdudgeon/fragment-network-merges.
"""

import itertools
import getpass
import pandas as pd
import numpy as np
from neo4j import GraphDatabase
from rdkit import Chem

try:
    password
except NameError:
    password = getpass.getpass()

driver = GraphDatabase.driver("bolt://localhost:7687", auth=("swills", password))

class Nodes():

    def __init__(self, fragments, names):
        self.fragments = fragments  # the list of smiles of the fragments
        self.names = names  # the list of names of the fragments

    def find_molecule_node(self, tx, smiles):
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

    def check_for_nodes(self):
        """
        Checks if nodes are present.
        Prints statement with how many of the fragments are present in the network.
        """
        number_nodes = 0
        number_fragments = len(self.fragments)
        with driver.session() as session:
            for fragment in self.fragments:
                mol = session.read_transaction(self.find_molecule_node, fragment)  # run neo4j query
                if mol:
                    number_nodes += 1
            print(f'{number_nodes} of {number_fragments} fragments present in network')

    def filter_for_nodes(self):
        """
        Checks if nodes are present and filters fragments for those that are in network.

        :param fragments: list of fragment smiles
        :type fragments: list

        :return: list of available fragments
        :rtype: list of smiles strings
        :return: list of available fragments (names)
        :rtype: list of strings
        """
        removed = 0
        with driver.session() as session:
            for i, fragment in enumerate(self.fragments):
                mol = session.read_transaction(self.find_molecule_node, fragment)  # run neo4j query
                if not mol:
                    self.fragments.pop(i)
                    self.names.pop(i)
                    removed += 1
        print(f'{removed} fragments removed from list. {len(self.fragments)} fragments remaining.')
        return self.fragments, self.names

class Merge():

    def __init__(self, fragments, names):
        self.fragments = fragments  # the list of smiles of the fragments
        self.names = names  # the list of names of the fragments
        self.fragment_pairs = []  # to store all permutations of fragments
        self.name_pairs = []  # to store all permutations of fragment names

    def add_required_synthons(self, labels, synthon):
        """
        Checks that the synthons have a single attachment point and single component.

        :param labels: set containing the synthons
        :type labels: set
        :param synthon: smiles string containing synthon with Xe atom
        :type synthon: string
        """
        if synthon.count('[Xe]') > 0:
            labels.add(synthon)

    def find_synthons(self, tx, smiles):
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
                self.add_required_synthons(labels, tokens[1])
                self.add_required_synthons(labels, tokens[4])
        return list(labels)

    def get_synthons(self, smiles):
        """
        Extract the synthons from the database.

        :param smiles: smiles of the fragment
        :type smiles: string

        :return: synthons
        :rtype: string
        """
        with driver.session() as session:
            synthons = session.read_transaction(self.find_synthons, smiles)
            print(f'Found {len(synthons)} synthons')
            return synthons

    def find_expansions(self, tx, smiles, synthon):
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

    def get_expansions(self, fragmentA, fragmentB):
        """
        Get the expansions of one fragment using a set of synthons from the other fragment.
        Returns a list of the expansions and a list of the synthons used to generate
        the expansions (this information is useful for later filtering).

        :param fragment: smiles string of the fragment to expand
        :type fragment: string
        :param synthons: list of synthons smiles
        :type synthons: list

        :return: list containing expansions
        :rtype: list of smiles
        :return: list containing synthons used
        :rtype: list of smiles
        """
        synthons = self.get_synthons(fragmentB)
        all_expansions = []
        all_synthons = []

        with driver.session() as session:
            number = 0
            total_expansions = 0
            expanded_synthons = 0
            for synthon in synthons:
                expansions = session.read_transaction(self.find_expansions, fragmentA, synthon)
                for smi in expansions:
                    all_expansions.append(smi)
                    all_synthons.append(synthon)
                print(f'Synthon {number}: found {len(expansions)} expansions')
                number += 1
                total_expansions += len(expansions)
                if expansions:  # record if the synthon led to expansions
                    expanded_synthons += 1
        print(f'{total_expansions} expansions from {expanded_synthons} out of {len(synthons)} synthons')
        print('')
        return all_expansions, all_synthons

    def get_combinations(self):
        """
        Enumerate all possible combinations of fragments for merging.

        :return: list of possible combinations of fragments
        :rtype: list of tuples
        :return: list of possible combinations of fragments (names)
        :rtype: list of tuples
        """
        self.fragment_pairs = list(itertools.permutations(self.fragments, 2))
        self.name_pairs = list(itertools.permutations(self.names, 2))
        return self.fragment_pairs, self.name_pairs

    def get_all_expansions(self):
        """
        Function executes the whole process, enumerating all possible fragment pairs,
        generating synthons for fragment B and generating all possible expansions.
        Returns a dataframe containing all the information, including the parent fragments, 
        the synthon used and the smiles of the merge.

        :return: dataframe containing all the expansions and related info
        :rtype: pandas dataframe
        """
        self.get_combinations()  # enumerate all fragment pairs
        col_names = ['Merge', 'Fragment A ID', 'Fragment A smiles', 'Fragment B ID', 'Fragment B smiles', 'Synthon', 'Smiles']
        df = pd.DataFrame()  # create empty df to store data
        for pair, name in zip(self.fragment_pairs, self.name_pairs):
            fragmentA, fragmentB = pair[0], pair[1]
            nameA, nameB = name[0], name[1]
            print(f'Expanding fragment A: {nameA} with synthons of fragment B: {nameB}')
            expansions, synthons = self.get_expansions(fragmentA, fragmentB)  # generate expansions
            # create dataframe to store results
            merge = nameA + '_' + nameB
            merges = [merge] * len(expansions)
            A_col = [fragmentA] * len(expansions)
            A_col_id = [nameA] * len(expansions)
            B_col = [fragmentB] * len(expansions)
            B_col_id = [nameB] * len(expansions)
            data = list(zip(merges, A_col_id, A_col, B_col_id, B_col, synthons, expansions))
            sub_df = pd.DataFrame(data, columns=col_names)
            df = df.append(sub_df, ignore_index=True)  # append df for each pair
        return df

# test data doesn't enumerate all combinations - use class to specify pairs
class MergePairs():

    def __init__(self, fragmentAs, fragmentBs, namesA, namesB):
        self.fragmentAs = fragmentAs  # the list of smiles of the fragments
        self.fragmentBs = fragmentBs  # the list of smiles of the fragments
        self.namesA = namesA  # the list of names of the fragments
        self.namesB = namesB  # the list of names of the fragments
        self.fragment_pairs = []  # to store all permutations of fragments
        self.name_pairs = []  # to store all permutations of fragment names

    def add_required_synthons(self, labels, synthon):
        """
        Checks that the synthons have a single attachment point and single component.

        :param labels: set containing the synthons
        :type labels: set
        :param synthon: smiles string containing synthon with Xe atom
        :type synthon: string
        """
        if synthon.count('[Xe]') > 0:
            labels.add(synthon)

    def find_synthons(self, tx, smiles):
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
                self.add_required_synthons(labels, tokens[1])
                self.add_required_synthons(labels, tokens[4])
        return list(labels)

    def get_synthons(self, smiles):
        """
        Extract the synthons from the database.

        :param smiles: smiles of the fragment
        :type smiles: string

        :return: synthons
        :rtype: string
        """
        with driver.session() as session:
            synthons = session.read_transaction(self.find_synthons, smiles)
            print(f'Found {len(synthons)} synthons')
            return synthons

    def find_expansions(self, tx, smiles, synthon):
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

    def get_expansions(self, fragmentA, fragmentB):
        """
        Get the expansions of one fragment using a set of synthons from the other fragment.
        Returns a list of the expansions and a list of the synthons used to generate
        the expansions (this information is useful for later filtering).

        :param fragment: smiles string of the fragment to expand
        :type fragment: string
        :param synthons: list of synthons smiles
        :type synthons: list

        :return: list containing expansions
        :rtype: list of smiles
        :return: list containing synthons used
        :rtype: list of smiles
        """
        synthons = self.get_synthons(fragmentB)
        all_expansions = []
        all_synthons = []

        with driver.session() as session:
            number = 0
            total_expansions = 0
            expanded_synthons = 0
            for synthon in synthons:
                expansions = session.read_transaction(self.find_expansions, fragmentA, synthon)
                for smi in expansions:
                    all_expansions.append(smi)
                    all_synthons.append(synthon)
                print(f'Synthon {number}: found {len(expansions)} expansions')
                number += 1
                total_expansions += len(expansions)
                if expansions:  # record if the synthon led to expansions
                    expanded_synthons += 1
        print(f'{total_expansions} expansions from {expanded_synthons} out of {len(synthons)} synthons')
        print('')
        return all_expansions, all_synthons

    def get_combinations(self):
        """
        Get fragment pairs and names as list of tuples.

        :return: list of possible combinations of fragments
        :rtype: list of tuples
        :return: list of possible combinations of fragments (names)
        :rtype: list of tuples
        """
        for fA, fB in zip(self.fragmentAs, self.fragmentBs):
            self.fragment_pairs.append((fA, fB))
        for nameA, nameB in zip(self.namesA, self.namesB):
            self.name_pairs.append((nameA, nameB))
        return self.fragment_pairs, self.name_pairs

    def get_all_expansions(self):
        """
        Function executes the whole process, enumerating all possible fragment pairs,
        generating synthons for fragment B and generating all possible expansions.
        Returns a dataframe containing all the information, including the parent fragments, 
        the synthon used and the smiles of the merge.

        :return: dataframe containing all the expansions and related info
        :rtype: pandas dataframe
        """
        self.get_combinations()  # enumerate all fragment pairs
        col_names = ['Merge', 'Fragment A ID', 'Fragment A smiles', 'Fragment B ID', 'Fragment B smiles', 'Synthon', 'Smiles']
        df = pd.DataFrame()  # create empty df to store data
        for pair, name in zip(self.fragment_pairs, self.name_pairs):
            fragmentA, fragmentB = pair[0], pair[1]
            nameA, nameB = name[0], name[1]
            print(f'Expanding fragment A: {nameA} with synthons of fragment B: {nameB}')
            expansions, synthons = self.get_expansions(fragmentA, fragmentB)  # generate expansions
            # create dataframe to store results
            merge = nameA + '_' + nameB
            merges = [merge] * len(expansions)
            A_col = [fragmentA] * len(expansions)
            A_col_id = [nameA] * len(expansions)
            B_col = [fragmentB] * len(expansions)
            B_col_id = [nameB] * len(expansions)
            data = list(zip(merges, A_col_id, A_col, B_col_id, B_col, synthons, expansions))
            sub_df = pd.DataFrame(data, columns=col_names)
            df = df.append(sub_df, ignore_index=True)  # append df for each pair
        return df
