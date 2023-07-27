"""
Used to generate merges between two fragments using the fragment network using port-forwarding.
Uses code from https://github.com/tdudgeon/fragment-network-merges.
"""

import getpass

from fragment_network_merges.merge.config_merge import config_merge
from fragment_network_merges.merge.find_merges_generic import (
    MergerFinder_generic,
    SearchSession_generic,
    add_required_synthons,
)
from neo4j import GraphDatabase


class Neo4jDriverWrapper(SearchSession_generic):
    def __init__(self, session):
        self.session = session

    def __enter__(self):
        """
        Starts the neo4j search session
        """
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        """
        Closes the neo4j search session
        """
        self.session.close()

    def _find_molecule_node(self, tx, smiles: str):
        """
        Finds node in the fragment network.

        :param tx: transaction in which query is run
        :type tx: neo4j transaction
        :param smiles: smiles of the fragment
        :type smiles: string

        :return: molecule node
        :rtype: node
        """
        for record in tx.run("MATCH (m:F2 {smiles: $smiles}) RETURN m", smiles=smiles):
            node = record["m"]
            return node

    def find_molecule_node(self, fragment: str):
        """
        Implements self._find_molecule_node()
        """
        return self.session.read_transaction(self._find_molecule_node, fragment)

    def _find_synthons(self, tx, smiles: str) -> list:
        """
        Query for all child fragments (recursive).
        Extract the label property of each edge and collect a set of SMILES that match our needs.

        :param tx: neo4j transaction
        :type tx: neo4j transaction
        :param smiles: smiles string of molecule to fragment
        :type smiles: str

        :returns: list of synthons
        :rtype: list
        """
        labels = set()
        for record in tx.run(  # if set to unlimited number of hops, the queries seem to stall
            "MATCH (fa:F2 {smiles: $smiles})-[e:FRAG*0..20]->(f:F2) RETURN e", smiles=smiles
        ):
            edges = record["e"]
            for edge in edges:
                s = edge["label"]
                tokens = s.split("|")
                add_required_synthons(labels, tokens[1])
                add_required_synthons(labels, tokens[4])
        return list(labels)

    def find_synthons(self, fragment: str):
        """
        Implements self._find_synthons()
        """
        return self.session.read_transaction(self._find_synthons, fragment)

    def _find_expansions(
        self,
        tx,
        smiles: str,
        synthon: str,
        number_hops: int,
        min_hac: int,
        max_hac: int,
        min_hac_fa: int,
    ) -> set:
        """
        Expand fragment 'A' using the synthons generated from fragment 'B' using a neo4j
        query. Query limited to compounds available from vendors a specified number of hops away,
        and filtered according to min/max heavy atom count.

        :param smiles: smiles of the fragment to expand
        :type smiles: str
        :param synthon: synthon of the generated synthon
        :type synthon: str
        :param number_hops: number of hops away from fragment A to look for merges
        :type number_hops: int
        :param min_hac: minimum number of heavy atoms of merges
        :type min_hac: int
        :param max_hac: maximum number of heavy atoms of merges
        :type max_hac: int
        :param min_hac_fa: minimum number of heavy atoms of the node before expansion
        :type min_hac_fa: int

        :return: expansions
        :rtype: set
        """
        synthon_no_xe = synthon.replace("([Xe])", "")
        synthon_no_xe = synthon_no_xe.replace("[Xe]", "")

        query = (
            "MATCH (fa:F2 {smiles: $smiles})"
            "-[:FRAG*0..%(number_hops)d]-(fb:F2)"
            "<-[e:FRAG]-(c:Mol) WHERE"
            " NOT fb.smiles = '%(synthon_no_xe)s'"  # check node before expansion is not equiv to synthon
            " AND fb.hac >= %(min_hac_fa)d"  # check num heavy atoms of node before expansion > 5
            " AND %(min_hac)d <= c.hac <= %(max_hac)d AND"  # check heavy atoms of final mol
            " (split(e.label, '|')[1] = $synthon OR split(e.label, '|')[4] = $synthon)"  # expanded with synthon
            " RETURN DISTINCT c"
            % {
                "number_hops": number_hops,
                "min_hac": min_hac,
                "max_hac": max_hac,
                "synthon_no_xe": synthon_no_xe,
                "min_hac_fa": min_hac_fa,
            }
        )
        expansions = set()
        for record in tx.run(query, smiles=smiles, synthon=synthon):
            node = record["c"]
            expansions.add(node["smiles"])
        return expansions

    def find_expansions(
        self,
        smiles: str,
        synthon: str,
        number_hops: int = config_merge.NUM_HOPS,
        min_hac: int = config_merge.MIN_HAC,
        max_hac: int = config_merge.MAX_HAC,
        min_hac_fa: int = config_merge.MIN_HAC_FA,
    ) -> set:
        """
        Implements self.find_expansions
        """
        return self.session.read_transaction(
            self._find_expansions,
            smiles,
            synthon,
            number_hops,
            min_hac,
            max_hac,
            min_hac_fa,
        )

    def _add_substituent(
        self,
        tx,
        smiles: str,
        substituent: str,
    ) -> set:

        query = (
                "MATCH (fa:F2 {smiles: $smiles})"
                "<-[e:FRAG]-(c:Mol) WHERE"
                " (split(e.label, '|')[1] = $substituent OR split(e.label, '|')[4] = $substituent)"
                " RETURN DISTINCT c"
        )
        expansions = set()
        for record in tx.run(query, smiles=smiles, substituent=substituent):
            node = record["c"]
            expansions.add(node["smiles"])
        return expansions

    def add_substituent(
        self,
        smiles: str,
        substituent: str,
    ) -> set:
        """
        Implements self.find_expansions
        """
        return self.session.read_transaction(
            self._add_substituent,
            smiles,
            substituent
        )


class MergerFinder_neo4j(MergerFinder_generic):
    def __init__(self, **kwargs):
        self._driver = None

    @property
    def driver(self):
        if self._driver is None:
            if not config_merge.NEO4J_PASS:
                password = getpass.getpass()
            else:
                password = config_merge.NEO4J_PASS
            self._driver = GraphDatabase.driver(
                config_merge.NEO4J_URI,# "bolt://localhost:7687",
                auth=(config_merge.NEO4J_USER, password)
            )
        return self._driver

    @driver.deleter
    def driver(self):
        self._driver.close()
        self._driver = None

    def getSearchSession(self):
        return Neo4jDriverWrapper(self.driver.session())
