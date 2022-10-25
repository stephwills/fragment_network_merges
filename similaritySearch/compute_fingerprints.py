import os

import numpy as np
from rdkit import Chem
from rdkit import DataStructs
from rdkit.Chem import rdMolDescriptors, AllChem
from rdkit.DataStructs import cDataStructs, SparseBitVect
from rdkit.Chem.Pharm2D import Generate

import similaritySearch.similaritySearchConfig as config

def get_fingerprint(smi):
    mol = Chem.MolFromSmiles(smi)
    if mol == None:
        return None
    return get_fingerprint_function()(mol)

def _morgan_fingerprint(mol):
    finprint = rdMolDescriptors.GetMorganFingerprintAsBitVect(
        mol, radius=config.FINGERPRINT_RADIUS, nBits=config.FINGERPRINT_NBITS, useChirality=config.FINGERPRINT_USE_CHIRALITY,
        useBondTypes=config.FINGERPRINT_USE_BOND_TYPES, useFeatures=config.FINGERPRINT_USE_FEATURES)
    return finprint

def _pharmacophore_fingerprint(mol):
    finprin = Generate.Gen2DFingerprint(mol, _get_sigFactory())
    onBits = [elem%config.FINGERPRINT_NBITS for elem in finprin.GetOnBits()]
    folded_fp = SparseBitVect(config.FINGERPRINT_NBITS)
    folded_fp.SetBitsFromList(onBits)
    return folded_fp

SIG_FACTORY=None
def _get_sigFactory():
    from rdkit.Chem.Pharm2D.SigFactory import SigFactory
    from rdkit import RDConfig
    global SIG_FACTORY
    if SIG_FACTORY is None:
        # featFactory = AllChem.BuildFeatureFactory(os.path.join(RDConfig.RDDataDir, 'BaseFeatures.fdef'))
        featFactory = AllChem.BuildFeatureFactory(os.path.join(os.path.dirname(__file__), 'data/pharmaFeatures.fdef'))

        SIG_FACTORY= SigFactory(featFactory,minPointCount=2,maxPointCount=3)
        SIG_FACTORY.SetBins([(0, 1), (1, 2), (2, 3), (3, 4), (4, 5), (5, 8), (8, 100)])
        SIG_FACTORY.Init()
        # assert SIG_FACTORY.GetSigSize() == config.FINGERPRINT_NBITS, f"Error, current ph4 fingerprint has {SIG_FACTORY.GetSigSize()} bits but config says {config.FINGERPRINT_NBITS}"
    return SIG_FACTORY

FINGERPRINT_FUN = None
def get_fingerprint_function():
    global FINGERPRINT_FUN
    if FINGERPRINT_FUN is None:
        if config.FINGERPRINT_TYPE == "morgan":
            FINGERPRINT_FUN = _morgan_fingerprint
        elif config.FINGERPRINT_TYPE == "pharmacophore":
            FINGERPRINT_FUN = _pharmacophore_fingerprint
        else:
            raise NotImplementedError("FINGERPRINT_TYPE is not valid. Only morgan and pharmacophore implemented")
    return FINGERPRINT_FUN

def computeFingerprintStr(smi):
    # print(smi)
    finPrint = get_fingerprint(smi)
    if finPrint is None:
        return b""
    finStr = cDataStructs.BitVectToBinaryText(finPrint)
    return finStr

#This is the one we use
def get_fingerPrint_as_npBool(smi):
    mol = Chem.MolFromSmiles(smi)
    if mol == None:
        return None
    finprin= rdMolDescriptors.GetMorganFingerprintAsBitVect(
        mol, radius=2, nBits=config.FINGERPRINT_NBITS, useChirality=0, useBondTypes=1, useFeatures=0)

    fp_num = np.zeros((0,), dtype=np.bool)
    DataStructs.ConvertToNumpyArray(finprin, fp_num)
    return fp_num



def decompressFingerprint_npStr(fpr_np):
    return np.unpackbits(np.frombuffer(fpr_np, dtype=np.uint8)).astype(bool)
