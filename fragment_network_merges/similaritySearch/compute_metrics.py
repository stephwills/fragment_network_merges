import numba
import numpy as np

from fragment_network_merges.similaritySearch import similaritySearchConfig as config


def jaccard_vectorized(x,y):
  intersections_count= x.astype(np.float32) @ y.T
  counts_x = np.sum(x, axis=-1)
  counts_y = np.sum(y, axis=-1)
  sums = counts_x.reshape(-1,1)+counts_y.reshape(-1,1).T
  union = sums - intersections_count
  return intersections_count / union


@numba.jit(**config.NUMBA_kwARGS)
def jaccard_numba(query_fp, db_fp):

    n11 = np.sum((query_fp == 1) & (db_fp == 1))
    n00 = np.sum((query_fp == 0) & (db_fp == 0))
    jac = n11 / (query_fp.shape[0] - n00)
    return jac


@numba.jit(**config.NUMBA_kwARGS)
def tversky_numba(query_fp, db_fp, alpha=0.3, beta=0.7):

    query_1 = (query_fp == 1)
    n11 = np.sum( query_1 & (db_fp == 1))
    n10 = np.sum( query_1 & (db_fp == 0))
    n01 = np.sum((query_fp == 0) & (db_fp == 1))
    res = n11 / (n11 + alpha*n10 +beta*n01)
    return res

@numba.jit(**config.NUMBA_kwARGS)
def fraction_of_query_on_bits(query_fp, db_fp):

    n11 = np.sum((query_fp == 1) & (db_fp == 1))
    score = n11 / np.sum(query_fp)
    return score

@numba.jit(**config.NUMBA_kwARGS)
def one_minus_fraction_of_query_on_bits(query_fp, db_fp):

    n11 = np.sum((query_fp == 1) & (db_fp == 1))
    score = 1 - n11 / np.sum(query_fp)
    return score

@numba.jit(**config.NUMBA_kwARGS)
def jaccard_weighted_numba(query_fp, db_fp, log_fp_frequency): #https://arxiv.org/pdf/1902.03402.pdf

    fp = np.stack((query_fp, db_fp), 0).astype(np.float64)
    fp = fp * log_fp_frequency
    # jac = np.sum(fp.min(axis=0)/fp.max(axis=0))
    jac = 0.
    for i in range(query_fp.shape[0]):     ### np.min(x, axis=0) is not supported, thus the loop
        numerator = min(fp[0,i], fp[1,i])
        denominator = max(fp[0,i], fp[1,i])
        if denominator > 0.:
            jac += numerator/denominator

    return jac

@numba.jit(**config.NUMBA_kwARGS)
def fp_bits_frequency(query_fp, db_fp, bits_logprob):
    positive_bits = (query_fp == 1) & (db_fp == 1)
    return np.sum(positive_bits / bits_logprob)


@numba.njit(cache=True, parallel=True, nogil=True)
def numba_logsumexp_stable(p):
    n, m = p.shape
    out = np.zeros(n)
    for i in numba.prange(n):
        p_max = np.max(p[i])
        res = 0
        for j in range(m):
            res += np.exp(p[i, j] - p_max)
        res = np.log(res) + p_max
        out[i] = res
    return out

@numba.jit(**config.NUMBA_kwARGS)
def fp_bits_frequency2D(query_fp, db_fp, bits_logprob2D):
    """
    Assuming only depth 1 in conditional dependencies: P(bi) = sum_i P(bi|bj)
    db_fp_frequency_weights := P(bi)
    db_fp_frequency2D_weights := P(bi,bj) = P(bi|bj)P(bj)
    P(bi|bj) := db_fp_frequency2D_weights / db_fp_frequency_weights
    P(bi) = ((positive_bits*positive_bits.T) * (db_fp_frequency2D_weights / db_fp_frequency_weights)).sum(0)
    """
    positive_bits = (query_fp == 1) & (db_fp == 1)
    postive_bits_mask = positive_bits.reshape(1,-1) & positive_bits.reshape(-1,1)
    postive_bits_mask = postive_bits_mask
    # with numba.objmode(logprob_bits='float64[:]'): #TODO: Implement logsumexp in numba to speed up code
        # logprob_bits = logsumexp(postive_bits_mask*bits_logprob2D, axis=0, return_sign=False)
    logprob_bits = numba_logsumexp_stable(postive_bits_mask * bits_logprob2D)
    return logprob_bits.sum()

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
  from fragment_network_merges.similaritySearch.compute_fingerprints import get_fingerprint
  fp1 = get_fingerprint(f1_smi)
  fp2 = get_fingerprint(f2_smi)
  fp_merge = get_fingerprint(merge_smi)

  print("smi1, smi2", DataStructs.FingerprintSimilarity(fp1, fp2))
  print("smi1, merge_smi", DataStructs.FingerprintSimilarity(fp1, fp_merge))
  print("smi2, merge_smi", DataStructs.FingerprintSimilarity(fp2, fp_merge))


def testTversky():
  from rdkit import DataStructs
  from fragment_network_merges.similaritySearch.compute_fingerprints import get_fingerprint

  f1_smi, f2_smi, merge_smi = _getTestInput()

  fp_merge = get_fingerprint(merge_smi)
  fp1 = get_fingerprint(f1_smi)
  fp2 = get_fingerprint(f2_smi)

  tversky_params = (0.3, 0.7)
  sim1 = DataStructs.FingerprintSimilarity(fp_merge, fp1,
                                           metric=lambda fp1, fp2: DataStructs.TverskySimilarity(fp1, fp2,
                                                                                                 *tversky_params))
  print("merge vs 1:", sim1)
  from fragment_network_merges.similaritySearch.compute_fingerprints import get_fingerPrint_as_npBool
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