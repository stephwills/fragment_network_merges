"""
Used to generate merges between two fragments using the fragment network.
Uses code from https://github.com/tdudgeon/fragment-network-merges.
"""

import os
import json
import itertools

from abc import ABC, abstractmethod
from typing import List, Tuple
from rdkit import Chem
from rdkit.Chem import rdShapeHelpers, Mol

from merge.config_merge import config_merge
from merge.preprocessing import get_mol
from filter.embedding_filter import add_coordinates, remove_xe


class SearchSession_generic(ABC):
    """
    Abstract class for searching the Fragment Network.
    """
    @abstractmethod
    def find_synthons(self, smiles:str) -> List[str]:
        """
        To extract synthons for a given fragment node.
        """
        raise NotImplementedError()

    @abstractmethod
    def find_molecule_node(self, fragment:str):
        """
        Retrieves a node for a given fragment from the database.
        """
        raise NotImplementedError()

    @abstractmethod
    def find_expansions( self, fragmentA:str, synthon:str) -> set:
        """
        Identifies expansions of a fragment that contain a synthon of interest.
        """
        raise NotImplementedError()


class MergerFinder_generic(ABC):
    """
    Abstract class containing functions for finding and filtering fragment nodes, filtering synthons and finding
    expansions
    """
    @abstractmethod
    def getSearchSession(self) -> SearchSession_generic:
        """
        Get neo4j search session.
        """
        raise NotImplementedError()

    def check_for_nodes(self, fragments:list) -> Tuple[int, int]:
        """
        Checks if nodes are present.
        Prints statement with how many of the fragments are present in the network.

        :param fragments: list of smiles of the fragments
        :type fragments: list

        :return: number of nodes present in network
        :rtype: int
        :return: number of fragments queried
        :rtype: int
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

    def filter_for_nodes(self, fragments:list, names:list) -> Tuple[list, list]:
        """
        Checks if nodes are present and filters fragments for those that are in network.

        :param fragments: list of fragment smiles
        :type fragments: list
        :param names: list of fragment names
        :type names: list

        :return: list of available fragments as smiles strings
        :rtype: list
        :return: list of available fragments (names, e.g. x0034)
        :rtype: list
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

    def get_synthons(self, smiles:str) -> list:
        """
        Extracts the synthons from the database for a given fragment.

        :param smiles: smiles of the fragment
        :type smiles: string

        :return: list of synthons strings
        :rtype: list
        """
        with self.getSearchSession() as session:
            synthons = session.find_synthons(smiles)
            print(f'Found {len(synthons)} synthons')
            return synthons

    def get_combinations(self, fragments:list, names:list) -> Tuple[list, list]:
        """
        Enumerate all possible combinations of fragments for merging.

        :param fragments: list of fragment smiles
        :type fragments: list
        :param names: list of fragment names
        :type names: list

        :return: list of possible combinations (as tuples) of fragments
        :rtype: list
        :return: list of possible combinations (as tuples) of fragments (names)
        :rtype: list
        """
        fragment_pairs = list(itertools.permutations(fragments, 2))
        name_pairs = list(itertools.permutations(names, 2))
        return fragment_pairs, name_pairs

    def carbons_check(self, synthons:list, num_carbons:int=config_merge.MIN_CARBONS) -> list:
        """
        Counts the number of carbons in the synthons. Keeps synthons that have at least three carbons.

        :param synthons: list of synthons as smiles to filter
        :type synthon: list
        :param num_carbons: minimum number of carbons for the synthon to be allowed to merge
        :type num_carbons: int

        :return: list of filtered synthons as smiles
        :rtype: list
        """
        filtered_synthons = []
        carbon = Chem.MolFromSmarts('[#6]')
        for synthon in synthons:
            synthon_mol = Chem.MolFromSmiles(synthon)
            carbon_count = len(synthon_mol.GetSubstructMatches(carbon))
            if carbon_count >= num_carbons:
                filtered_synthons.append(synthon)

        return filtered_synthons

    def substructure_check(self, synthons:list, fragmentA:Mol, fragmentB:Mol, synthon_overlap:float=config_merge.SYNTH_OVERLAP)\
            -> list:
        """
        Checks if the synthon is already present in fragment A and in an overlapping position.
        These synthons result in elaborations of fragment A rather than a merge.

        :param synthons: list of synthons as smiles to filter
        :type synthon: list
        :param fragment A: fragment undergoing expansions
        :type fragment A: RDKit molecule
        :param fragment B: fragment to generate synthons
        :type fragment B: RDKit molecule
        :param synthon_overlap: the maximum amount of overlap between synthons to be classified as not a merge
        :type synthon_overlap: float

        :return: list of filtered synthons as smiles
        :rtype: list
        """
        filtered_synthons = []
        for synthon in synthons:
            # add coordinates from the fragments to the synthon structure
            synthon_A = Chem.MolFromSmiles(synthon)
            synthon_A = remove_xe(synthon_A)
            synthon_B = Chem.Mol(synthon_A)

            try:
                fA_match = add_coordinates(fragmentA, synthon_A)
                fB_match = add_coordinates(fragmentB, synthon_B)

                distance = rdShapeHelpers.ShapeProtrudeDist(fA_match, fB_match)
                protrusion_threshold = 1 - synthon_overlap
                if protrusion_threshold >= distance:
                    filtered_synthons.append(synthon)

            except:
                filtered_synthons.append(synthon)

        return filtered_synthons

    def get_expansions(self, fragments:tuple, names:tuple, target:str, output_dir:str) -> dict:
        """
        Function executes the whole process, generating synthons for fragment B and using them to
        generate expansions of fragment A. Returns a dictionary containing all the synthons as keys,
        and the smiles of the merges they were used to generate as values.
        The dictionary is then saved as a json file.

        :param fragments: tuple of fragment smiles for merging
        :type fragments: tuple
        :param names: tuple of fragment names for merging
        :type names: tuple
        :param target: name of protein target (e.g. 'nsp13')
        :type target: str
        :param output_dir: name of directory to save output files
        :type output_dir: str

        :return: dictionary of merges
        :rtype: dict
        """
        # get fragments A and B from the tuple
        fragmentA, fragmentB = fragments[0], fragments[1]
        nameA, nameB = names[0], names[1]
        molA, molB = get_mol(target, nameA), get_mol(target, nameB)
        print(f'Expanding fragment A: {nameA} with synthons of fragment B: {nameB}')

        # check if file already exists containing merges
        if output_dir is not None:
            filename = nameA + '_' + nameB + '.json'  # name json file using fragment names
            if not os.path.exists(output_dir):  # if output dir doesn't exist, create dir
                os.makedirs(output_dir)
            filepath = os.path.join(output_dir, filename)
            if os.path.isfile(filepath):
                with open(filepath) as f:
                    return json.load(f)  # if file exists, returns the existing merge dict

        # generate the synthons from fragment B
        synthons = self.get_synthons(fragmentB)

        # filter synthons
        synthons = self.carbons_check(synthons)  # filter by number of carbons
        synthons = self.substructure_check(synthons, molA, molB)  # filter by whether in fragment A
        print(f'{len(synthons)} synthons remaining after filtering')

        # run database query and expansion
        all_expansions = {}
        with self.getSearchSession() as session:
            number = 0  # keep track of synthon number
            total_expansions = 0
            expanded_synthons = 0
            for synthon in synthons:  # expand fragment A using each synthon
                print(f'Running synthon {number}')
                expansions = session.find_expansions(fragmentA, synthon)
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


def add_required_synthons(labels:set, synthon:str):
    """
    Checks that the synthon has a single attachment point and single component.
    If yes, the synthon is added to the set of labels.

    :param labels: set containing the synthons
    :type labels: set
    :param synthon: smiles string containing synthon with Xe atom
    :type synthon: string
    """
    if synthon.count('[Xe]') == 1 and synthon.count(".") == 0:
        labels.add(synthon)
