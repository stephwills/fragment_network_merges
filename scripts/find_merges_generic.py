"""
Used to generate merges between two fragments using the fragment network.
Uses code from https://github.com/tdudgeon/fragment-network-merges.
"""

import os
import itertools
import getpass
import json
from abc import ABC, abstractmethod

from rdkit import Chem
from rdkit.Chem import rdShapeHelpers

from scripts.preprocessing import get_mol
from scripts.embedding_filter import add_coordinates, remove_xe



class MergerFinder_generic(ABC):


    @abstractmethod
    def getSearchSession(self):
        raise NotImplementedError()
    
    # functions for checking the nodes and filtering for fragments that exist as nodes
    def check_for_nodes(self, fragments):
        """
        Checks if nodes are present.
        Prints statement with how many of the fragments are present in the network.

        :param fragments: list of smiles of the fragments
        :type fragments: list of strings
        """
        number_nodes = 0
        number_fragments = len(fragments)
        with self.getSearchSession() as session:
            for fragment in fragments:
                mol = session.find_molecule_node( fragment)  # run neo4j query

                if mol:
                    number_nodes += 1
            print(f'{number_nodes} of {number_fragments} fragments present in network')
        return number_nodes, number_fragments

    def filter_for_nodes(self, fragments, names):
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
        with self.getSearchSession() as session:
            for i, fragment in enumerate(fragments):
                mol = session.find_molecule_node(fragment)  # run neo4j query
                if not mol:
                    fragments.pop(i)
                    names.pop(i)
                    removed += 1
        print(f'{removed} fragments removed from list. {len(fragments)} fragments remaining.')
        return fragments, names

    def get_synthons(self, smiles):
        """
        Extract the synthons from the database.

        :param smiles: smiles of the fragment
        :type smiles: string

        :return: synthons
        :rtype: string
        """
        with self.getSearchSession() as session:
            synthons = session.find_synthons(smiles)
            print(f'Found {len(synthons)} synthons')
            return synthons

    def get_expansions(self, fragments, names, target, output_dir):
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
        if output_dir is not None:
            filename = nameA + '_' + nameB + '.json'
            if not os.path.exists(output_dir):
                os.makedirs(output_dir)
            filepath = os.path.join(output_dir, filename)
            if os.path.isfile(filepath):
                with open(filepath) as f:
                    return json.load(f)
        molA, molB = get_mol(target, nameA), get_mol(target, nameB)
        print(f'Expanding fragment A: {nameA} with synthons of fragment B: {nameB}')

        # generate the synthons from fragment B
        unfiltered_synthons = self.get_synthons(fragmentB)

        # filter synthons for those with <3 carbons
        synthons_3c = [self.filter_synthons(syn) for syn in unfiltered_synthons]
        # remove None values from list
        synthons_3c = list(filter(None, synthons_3c))

        # filter synthons for those already in fragment A
        synthons = [self.substructure_check(syn, molA, molB) for syn in synthons_3c]
        # remove None values from list
        synthons = list(filter(None, synthons))
        print(f'{len(synthons)} synthons remaining after filtering')

        # create empty dictionary to store results
        all_expansions = {}
        with self.getSearchSession() as session:
            number = 0
            total_expansions = 0
            expanded_synthons = 0
            for synthon in synthons:  # expand fragment A using each synthon
                print(f'Running synthon {number}')
                expansions = session.find_expansions( fragmentA, synthon)
                all_expansions[synthon] = list(expansions)  # store in dictionary with the synthon as key
                print(f'Synthon {number}: found {len(expansions)} expansions')
                number += 1
                total_expansions += len(expansions)
                if expansions:  # record if the synthon led to expansions
                    expanded_synthons += 1
        print(f'{total_expansions} expansions from {expanded_synthons} out of {len(synthons)} synthons')
        print('\n')

        # save as json file
        if output_dir is not None:
            with open(filepath, 'w') as f:
                json.dump(all_expansions, f)

        return all_expansions

    def get_combinations(self, fragments, names):
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



    def filter_synthons(self, synthon):
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

    def substructure_check(self, synthon, fragmentA, fragmentB):
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
            if distance >= 0.5:
                return synthon
            else:
                return None

        except:
            return synthon

def add_required_synthons( labels, synthon):
    """
    Checks that the synthons have a single attachment point and single component.

    :param labels: set containing the synthons
    :type labels: set
    :param synthon: smiles string containing synthon with Xe atom
    :type synthon: string
    """
    if synthon.count('[Xe]') == 1 and synthon.count(".") == 0:
        labels.add(synthon)