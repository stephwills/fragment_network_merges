"""
Used to generate merges between two fragments using the fragment network using RestAPI.
Uses code from https://github.com/tdudgeon/fragment-network-merges.
"""

import time
import requests
import urllib.parse

from datetime import datetime, MINYEAR

from merge.config_merge import config_merge
from merge.find_merges_generic import MergerFinder_generic, add_required_synthons, SearchSession_generic
from merge.utils import Singleton


class RestAPI_wrapper(SearchSession_generic):
    __metaclass__ = Singleton
    FRAGNET_SERVER = "https://fragnet-search.xchem-dev.diamond.ac.uk/fragnet-search-api/fragnet-search/rest/v2/search/"
    TOKEN_SERVER = "https://squonk.it/auth/realms/squonk/protocol/openid-connect/token"

    def __init__(self, num_retry=10, sleep_time_retry=5, token_validity_secs=3000):
        self.num_retry = num_retry
        self.sleep_time_retry = sleep_time_retry
        self.token_validity_secs = token_validity_secs
        self._token = None
        self._last_token_time = datetime(MINYEAR, 1, 1)  # Arbitrary old time

    @property
    def token(self):

        if (datetime.now() - self._last_token_time).total_seconds() > self.token_validity_secs:
            print("Generating new token for restAPI")
            token = None
            for i in range(self.num_retry):
                r = requests.post(self.TOKEN_SERVER, data=dict(grant_type="password", client_id="fragnet-search-ui",
                                                               username=config_merge.RESTAPI_USER,
                                                               password=config_merge.RESTAPI_PASS))
                if r.ok:
                    token = r.json()["access_token"]
                else:
                    time.sleep(self.sleep_time_retry)
            if token is None:
                print(r)
                raise ConnectionError("Error, authentication token couldn't be generated")
            self._last_token_time = datetime.now()
            self._token = token
        return self._token

    def find_molecule_node(self, fragment:str):
        """
        Finds node in the fragment network.

        :param fragment: smiles of the fragment
        :type fragment: str

        :return: molecule node
        :rtype: node
        """
        fragment = urllib.parse.quote_plus(fragment)
        result = None
        for i in range(self.num_retry):
            r = requests.get(f"{self.FRAGNET_SERVER}/molecule/{fragment}",
                             headers={"Authorization": f"bearer {self.token}"})
            if r.ok:
                result = r.json()
            else:
                time.sleep(self.sleep_time_retry)
            if result is None:
                print(r)
                raise ConnectionError("Error, find_molecule_node %s failed" % fragment)
        return result

    def _find_oneMol_synthons(self, smi):
        """
        Retrieves synthons from the Fragment Network.
        """
        smi = urllib.parse.quote_plus(smi)
        result = None
        for i in range(self.num_retry):
            r = requests.get(f"{self.FRAGNET_SERVER}/fragments/{smi}",
                             headers={"Authorization": f"bearer {self.token}"})
            if r.ok:
                result = r.json()
            else:
                time.sleep(self.sleep_time_retry)
        if result is None:
            print(r)
            raise ConnectionError("Error, _find_oneMol_synthons %s failed" % smi)
        return result

    def find_synthons(self, fragment:str) -> list:
        """
        Query for all child fragments (recursive).
        Extract the label property of each edge and collect a set of SMILES that match our needs.

        :param fragment: smiles string of molecule to fragment
        :type fragment: str

        :returns: list of synthons
        :rtype: list
        """
        labels = set()
        allSynthons = self._find_oneMol_synthons(fragment)
        for synthon in allSynthons:
            add_required_synthons(labels, synthon)
        return list(labels)

    def find_expansions(self, smiles: str, synthon: str, number_hops: int = config_merge.NUM_HOPS,
                        hacMin: int = config_merge.MIN_HAC, hacMax: int = config_merge.MAX_HAC) -> set:
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
        :param hacMin: minimum number of heavy atoms of merges
        :type hacMin: int
        :param hacMax: maximum number of heavy atoms of merges
        :type hacMax: int

        :return: expansions
        :rtype: set
        """

        smiles = urllib.parse.quote_plus(smiles)
        synthon = urllib.parse.quote_plus(synthon)

        result = None
        for i in range(self.num_retry):
            r = requests.get(f"{self.FRAGNET_SERVER}/synthon-expand/{smiles}?synthon={synthon}&"
                             f"hops={number_hops}&hacMin={hacMin}&hacMax={hacMax}",
                             headers={"Authorization": f"bearer {self.token}"})
            if r.ok:
                result = r.json()
            else:
                time.sleep(self.sleep_time_retry)
        if result is None:
            print(r)
            raise ConnectionError("Error, _find_expansions (%s - %s) failed" % (smi, synthon))
        expansions = set([example["smiles"] for example in result])
        return expansions

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        pass


class MergerFinder_restAPI(MergerFinder_generic):

    def __init__(self, **kwargs):
        self._driver = None

    def getSearchSession(self):
        return RestAPI_wrapper()
