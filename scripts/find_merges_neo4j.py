"""
Used to generate merges between two fragments using the fragment network.
Uses code from https://github.com/tdudgeon/fragment-network-merges.
"""

import os
import getpass

from neo4j import GraphDatabase

from scripts.config import config
from scripts.find_merges_generic import MergerFinder_generic, add_required_synthons, SearchSession_generic


class Neo4jDriverWrapper(SearchSession_generic):

    def __init__(self, session):
        self.session = session

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.session.close()

    def _find_molecule_node(self, tx, smiles):
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

    def find_molecule_node(self, fragment):
        return self.session.read_transaction(self._find_molecule_node, fragment)

    def _find_synthons(self, tx, smiles):
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

    def find_synthons(self, fragment):
        return self.session.read_transaction(self._find_synthons, fragment)


    def _find_expansions(self, tx, smiles, synthon, number_hops=2):
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
                    "-[:FRAG*0..%(number_hops)d]-(:F2)"  # to increase the number of hops, change '0..2' to '0..3'
                    "<-[e:FRAG]-(c:Mol) WHERE"
                    " c.hac > 15 AND"
                    " (split(e.label, '|')[1] = $synthon OR split(e.label, '|')[4] = $synthon)"
                    " RETURN DISTINCT c"%{"number_hops":number_hops})
        expansions = set()
        for record in tx.run(query, smiles=smiles, synthon=synthon):
            node = record['c']
            expansions.add(node['smiles'])
        return expansions

    def find_expansions(self, smiles, synthon, number_hops=2):
        return self.session.read_transaction(self._find_expansions, smiles, synthon, number_hops)


class MergerFinder_neo4j(MergerFinder_generic):

    def __init__(self, **kwargs):
        self._driver = None

    @property
    def driver(self):
        if self._driver is None:
            if not config.NEO4J_PASS:
                password = getpass.getpass()
            else:
                password = config.NEO4J_PASS
            self._driver = GraphDatabase.driver("bolt://localhost:7687",
                                          auth=(config.NEO4J_USER, password))
        return self._driver

    @driver.deleter
    def driver(self):
        self._driver.close()
        self._driver = None

    def getSearchSession(self):
        return Neo4jDriverWrapper(self.driver.session())


if __name__ == "__main__":
    from rdkit import Chem
    smi = Chem.MolToSmiles(Chem.MolFromSmiles("CC=1C=CC(CS(=O)(=O)N)=CC1F"))
    print(smi)
    merger = MergerFinder_neo4j()
    results = merger.check_for_nodes([smi])
    print(results)

    results = merger.get_synthons(smi)
    print(results)

    smi2 = Chem.MolToSmiles(Chem.MolFromSmiles("CS(=O)(=O)NCCC=1C=CC=CC1"))
    results = merger.get_expansions([smi, smi2], ["x0034_0B", "x0176_0B"], "nsp13", None)
    print(results)
