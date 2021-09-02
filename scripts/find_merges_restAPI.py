"""
Used to generate merges between two fragments using the fragment network.
Uses code from https://github.com/tdudgeon/fragment-network-merges.
"""

import os
import getpass
import time
import urllib.parse
from datetime import datetime, MINYEAR

import requests

from scripts.config import config
from scripts.find_merges_generic import MergerFinder_generic, add_required_synthons
from scripts.utils import Singleton

class RestAPI_wrapper(metaclass=Singleton):

    FRAGNET_SERVER = "https://fragnet-search.xchem-dev.diamond.ac.uk/fragnet-search-api/fragnet-search/rest/v2/search/"
    TOKEN_SERVER = "https://squonk.it/auth/realms/squonk/protocol/openid-connect/token"

    def __init__(self, num_retry=10, sleep_time_retry=5, token_validity_secs=3000):
        self.num_retry = num_retry
        self.sleep_time_retry = sleep_time_retry
        self.token_validity_secs = token_validity_secs
        self._token = None
        self._last_token_time = datetime(MINYEAR,1,1) #Arbitrary old time

    @property
    def token(self):

        if (datetime.now() - self._last_token_time).total_seconds() > self.token_validity_secs:
            print("Generating new token for restAPI")
            token = None
            for i in range(self.num_retry):
                r = requests.post(self.TOKEN_SERVER, data=dict(grant_type="password", client_id="fragnet-search-ui",
                                                               username=config.RESTAPI_USER, password=config.RESTAPI_PASS))
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

    def find_molecule_node(self, fragment):

        fragment = urllib.parse.quote_plus( fragment )
        result = None
        for i in range(self.num_retry):
            r = requests.get(f"{self.FRAGNET_SERVER}/molecule/{fragment}", headers={"Authorization": f"bearer {self.token}"})
            if r.ok:
                result = r.json()
            else:
                time.sleep(self.sleep_time_retry)
            if result is None:
                print(r)
                raise ConnectionError("Error, find_molecule_node %s failed"%fragment)
        return result

    def _find_oneMol_synthons(self, smi):

        smi = urllib.parse.quote_plus( smi )
        result = None
        for i in range(self.num_retry):
            r =  requests.get(f"{self.FRAGNET_SERVER}/fragments/{smi}", headers={"Authorization": f"bearer {self.token}"})
            if r.ok:
                result = r.json()
            else:
                time.sleep(self.sleep_time_retry)
        if result is None:
            print(r)
            raise ConnectionError("Error, _find_oneMol_synthons %s failed"%smi)
        return result

    def find_synthons(self, fragment):
        labels = set()
        allSynthons =  self._find_oneMol_synthons(fragment)
        for synthon in allSynthons:
            add_required_synthons(labels, synthon)
        return list(labels)

    def find_expansions(self, smiles, synthon, number_hops=2, hacMin=3, hacMax=30):
        """
        Expand fragment 'A' using one synthon generated from fragment 'B' using a neo4j
        query. Query limited to compounds available from vendors, with HAC > 15
        and a maximum of 2 hops away.

        :param smiles: smiles of the fragment to expand
        :type smiles: string
        :param synthon: synthon of the generated synthon
        :type synthon: string

        :param hacMin: Minimun heavy atom count filter
        :type hacMin: int
        :param hacMax: Maximun heavy atom count filter
        :type hacMax: int

        :return: expansions
        :rtype: set
        """

        smiles = urllib.parse.quote_plus( smiles )
        synthon = urllib.parse.quote_plus( synthon )

        result = None
        for i in range(self.num_retry):
            r =  requests.get(f"{self.FRAGNET_SERVER}/synthon-expand/{smiles}?synthon={synthon}&"
                              f"hops={number_hops}&hacMin={hacMin}&hacMax={hacMax}",
                              headers={"Authorization": f"bearer {self.token}"})
            if r.ok:
                result = r.json()
            else:
                time.sleep(self.sleep_time_retry)
        if result is None:
            print(r)
            raise ConnectionError("Error, _find_expansions (%s - %s) failed"%( smi, synthon))
        expansions = set([ example["smiles"] for example in result ])
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


if __name__ == "__main__":
    from rdkit import Chem
    smi = Chem.MolToSmiles(Chem.MolFromSmiles("CC=1C=CC(CS(=O)(=O)N)=CC1F"))
    print(smi)
    merger = MergerFinder_restAPI()
    results = merger.check_for_nodes([smi])
    print(results)

    results = merger.get_synthons(smi)
    print(results)

    smi2 = Chem.MolToSmiles(Chem.MolFromSmiles("CS(=O)(=O)NCCC=1C=CC=CC1"))
    results = merger.get_expansions([smi, smi2], ["x0034_0B", "x0176_0B"], "nsp13", None)
    print(results)
