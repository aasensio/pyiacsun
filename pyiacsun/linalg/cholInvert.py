__all__ = ['cholInvert']
import numpy as np

def cholInvert(A):
    """
    Return the inverse of a symmetric matrix using the Cholesky decomposition. The log-determinant is
    also returned
    
    Args:
        A : (N,N) matrix
    
    Returns:
        AInv: matrix inverse
        logDeterminant: logarithm of the determinant of the matrix 
    """
    L = np.linalg.cholesky(A)
    LInv = np.linalg.inv(L)
    AInv = np.dot(LInv.T, LInv)
    logDeterminant = -2.0 * np.sum(np.log(np.diag(LInv)))   # Why the minus sign?
    return AInv, logDeterminant
