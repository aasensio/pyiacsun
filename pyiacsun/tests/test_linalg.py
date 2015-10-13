from __future__ import print_function

import numpy as np
import numpy.testing as npt
import pyiacsun

def test_fsvd():
    A = np.array([[1,2,3,4],[3,5,6,7],[2,8,3,1]])

    U2, S2, V2 = np.linalg.svd(A)
    U, S, V = pyiacsun.linalg.fsvd(A, 9, 3, usePowerMethod=True)

    npt.assert_equal(S, S2)

def test_cholinvert():
    A = np.array([[4,12,-16],[12,37,-43],[-16,-43,98]])

    AInv = np.linalg.inv(A)
    logDet = np.log(np.linalg.det(A))
    AInv2, logDet2 = pyiacsun.linalg.cholInvert(A)
    npt.assert_allclose(AInv, AInv2)
    npt.assert_allclose(logDet, logDet2)

def test_svht():
    rank = 3
    m = 80
    n = 80
    sigma = 0.02
    X1 = np.random.randn(m,rank)
    X2 = np.random.randn(rank,n)

    Y = np.dot(X1, X2) + sigma * np.random.randn(m,n)

    thrKnownNoise, thrUnknownNoise, noiseEstimation, wSVD = pyiacsun.linalg.optimalSVHT(Y)

    #print("Known sigma : {0} - rank={1}".format(thrKnownNoise * sigma, rankKnownNoise))
    #print("Unknown sigma : {0} - rank={1} - sigmaNoise={2}".format(thrUnknownNoise, rankUnknownNoise, noiseEstimation))
