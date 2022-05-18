import numba
import numpy as np


def jaccard_vectorized(x,y):
  intersections_count= x.astype(np.float32) @ y.T
  counts_x = np.sum(x, axis=-1)
  counts_y = np.sum(y, axis=-1)
  sums = counts_x.reshape(-1,1)+counts_y.reshape(-1,1).T
  union = sums - intersections_count
  return intersections_count / union


@numba.jit( nopython=True, cache=True)
def jaccard_numba(query_fp, db_fp):

    n11 = np.sum((query_fp == 1) & (db_fp == 1))
    n00 = np.sum((query_fp == 0) & (db_fp == 0))
    jac = n11 / (query_fp.shape[0] - n00)
    return jac


@numba.jit( nopython=True, cache=True)
def tversky_numba(query_fp, db_fp, alpha=0.3, beta=0.7):

    query_1 = (query_fp == 1)
    n11 = np.sum( query_1 & (db_fp == 1))
    n10 = np.sum( query_1 & (db_fp == 0))
    n01 = np.sum((query_fp == 0) & (db_fp == 1))
    res = n11 / (n11 + alpha*n10 +beta*n01)
    return res



def _getTestInput():
  # merge_smi = 'Cc1ccc([C@@](N)(O)OOCc2ccccc2)cc1F'
  # f1_smi = 'CC=1C=CC(CS(=O)(=O)N)=CC1'
  # f2_smi = 'OC=1C=CC(NC(=O)CCC=2C=CC=CC2)=CC1'

  # merge_smi = 'CN(C1CCCCCCC1)S(C)(=O)=O'
  # f1_smi = 'CN(C1CCCCCC1)S(=O)(=O)C'
  # f2_smi = 'NC=1C=CC(=CC1)S(=O)(=O)NC=2C=CC=CC2'


  merge_smi = 'CN(C1CCCCCCC1)S(C)(=O)=O'
  f1_smi =    'CN(C1CCCCCC1)S(=O)(=O)C'
  f2_smi =    'NC=1C=CC(=CC1)S(=O)(=O)NC=2C=CC=CC2'

  f1_smi, f2_smi = f2_smi, f1_smi

  return  f1_smi, f2_smi, merge_smi

def testTanimoto():
  from rdkit import DataStructs

  f1_smi, f2_smi, merge_smi = _getTestInput()
  from similaritySearch.compute_fingerprints import get_fingerPrint
  fp1 = get_fingerPrint(f1_smi)
  fp2 = get_fingerPrint(f2_smi)
  fp_merge = get_fingerPrint(merge_smi)

  print("smi1, smi2", DataStructs.FingerprintSimilarity(fp1, fp2))
  print("smi1, merge_smi", DataStructs.FingerprintSimilarity(fp1, fp_merge))
  print("smi2, merge_smi", DataStructs.FingerprintSimilarity(fp2, fp_merge))


def testTversky():
  from rdkit import DataStructs
  from similaritySearch.compute_fingerprints import get_fingerPrint

  f1_smi, f2_smi, merge_smi = _getTestInput()

  fp_merge = get_fingerPrint(merge_smi)
  fp1 = get_fingerPrint(f1_smi)
  fp2 = get_fingerPrint(f2_smi)

  tversky_params = (0.3, 0.7)
  sim1 = DataStructs.FingerprintSimilarity(fp_merge, fp1,
                                           metric=lambda fp1, fp2: DataStructs.TverskySimilarity(fp1, fp2,
                                                                                                 *tversky_params))
  print("merge vs 1:", sim1)
  from similaritySearch.compute_fingerprints import get_fingerPrint_as_npBool
  print( tversky_numba(get_fingerPrint_as_npBool(merge_smi), get_fingerPrint_as_npBool(f1_smi), *tversky_params) )

  sim2 = DataStructs.FingerprintSimilarity(fp_merge, fp2,
                                           metric=lambda fp1, fp2: DataStructs.TverskySimilarity(fp1, fp2,
                                                                                                 *tversky_params))
  print("merge vs 2:", sim2)
  print( tversky_numba(get_fingerPrint_as_npBool(merge_smi), get_fingerPrint_as_npBool(f2_smi), *tversky_params) )

if __name__ == "__main__":
  print("Tanimoto")
  testTanimoto()
  print("------------------------------")
  print("Tversky")
  testTversky()