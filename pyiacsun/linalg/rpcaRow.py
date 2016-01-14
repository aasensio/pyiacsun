"""Summary"""
from __future__ import print_function
__all__ = ['rpcaRow']

import numpy as np
from ipdb import set_trace as stop

def _proxLowRank(X, mu):
    """Summary
    
    Args:
        X (TYPE): Description
        mu (TYPE): Description
    
    Returns:
        TYPE: Description
    """
    U, w, V = np.linalg.svd(X, full_matrices=True)

    ind = np.where(np.abs(w) < mu)[0]
    w[ind] = 0.0
    ind = np.where(np.abs(w) > mu)[0]
    w[ind] = w[ind] - mu * np.sign(w[ind])

    n = V.shape[0]
    S = np.zeros((U.shape[1], V.shape[0]))
    S[:n,:n] = np.diag(w)

    return U.dot(S).dot(V)

def _proxSparse(X, mu):
    """Summary
    
    Args:
        X (TYPE): Description
        mu (TYPE): Description
    
    Returns:
        TYPE: Description
    """
    nRow, nCol = X.shape
    out = np.copy(X)

    for i in range(nRow):
        colNorm = np.linalg.norm(out[i,:], 2)
        if (colNorm < mu):
            out[i,:] = np.zeros(nCol)
        else:
            out[i,:] = out[i,:] - mu / colNorm * out[i,:]

    return out

def rpcaRow(X, lamb, Omega=None, tolerance=1e-5):
    """Compute the decomposition of the matrix X on the sum of two components
    
    X = L + S
    
    where L is low-rank and S is row-sparse. It is based on the paper "Robust PCA via Outlier Pursuit" (Xu et al. 2012)
    http://guppy.mpe.nus.edu.sg/~mpexuh/papers/OutlierPursuit-TIT.pdf
    http://guppy.mpe.nus.edu.sg/~mpexuh/publication.html

    Example:

        nl = 50
        nobs = 200
        rank = 3
        ncorr = 12

        x = np.random.randn(nobs,rank)
        y = np.random.randn(rank,nl)

        LObs = x.dot(y)
        CObs = np.zeros_like(LObs)

        ind = np.random.permutation(nobs)
        for i in range(ncorr):
            CObs[ind[i],:] = np.random.randn(nl)

        MObs = LObs + CObs #+ 1e-4 * np.random.randn(nobs, nl)

        L, S = rpcaRow(MObs, 0.35)
    
    Args:
        X (float): nxm matrix
        lamb (float): regularization parameter
        Omega (float, optional): matrix of active elements of the matrix. If not passed, all elements are assumed to be observed
        tolerance (float, optional): tolerance
    
    Returns:
        float: L (the low-rank component) and S (the row-sparse component)
    """

    if (Omega == None):
        Omega = np.ones_like(X)

    L0 = np.zeros_like(X)
    C0 = np.zeros_like(X)
    L1 = np.zeros_like(X)
    C1 = np.zeros_like(X)

    t0 = 1.0
    t1 = 1.0

    delta = 0.00001
    mu = 0.5 * np.linalg.norm(X)
    muBar = delta * mu
    eta = 0.9
    tol = tolerance * np.linalg.norm(X, 'fro')

    norm = 1e10

    loop = 0

    while (norm > tol**2):

        YL = L0 + (t1 - 1) / t0 * (L0 - L1)
        YC = C0 + (t1 - 1) / t0 * (C0 - C1)

        MDiff = (YL + YC - X) * Omega

        GL = YL - 0.5 * MDiff    
        LNew = _proxLowRank(GL, 0.5 * mu)

        GC = YC - 0.5 * MDiff    
        CNew = _proxSparse(GC, 0.5 * mu * lamb)

        t1 = np.copy(t0)
        t0 = 0.5 * (1.0 + np.sqrt(4.0 * t0**2 + 1))

        L1 = np.copy(L0)
        L0 = np.copy(LNew)

        C1 = np.copy(C0)
        C0 = np.copy(CNew)

        mu = np.max([eta * mu, muBar])

        SL = 2.0 * (YL - LNew) + (LNew + CNew - YL - YC)
        SC = 2.0 * (YC - CNew) + (LNew + CNew - YL - YC)

        norm = np.linalg.norm(SL, 'fro')**2 + np.linalg.norm(SC, 'fro')**2

        print("it: {0} - norm={1}".format(loop, norm))

        loop += 1

    return LNew, CNew


# Example

# nl = 50
# nobs = 200
# rank = 3
# ncorr = 12

# x = np.random.randn(nobs,rank)
# y = np.random.randn(rank,nl)

# LObs = x.dot(y)
# CObs = np.zeros_like(LObs)

# ind = np.random.permutation(nobs)
# for i in range(ncorr):
#     CObs[ind[i],:] = np.random.randn(nl)

#     MObs = LObs + CObs #+ 1e-4 * np.random.randn(nobs, nl)
# pl.close('all')
# f, ax = pl.subplots(ncols=3, nrows=2, figsize=(8,18))

# ax[0,0].imshow(LObs)
# ax[1,0].imshow(LNew)

# ax[0,1].imshow(CObs)
# ax[1,1].imshow(CNew)

# ax[0,2].imshow(MObs)
# ax[1,2].imshow(LNew + CNew)