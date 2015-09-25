__all__ = ['cholInvert']
import numpy as np
###########################################
# Return the inverse of a matrix and its logDeterminant
###########################################
def cholInvert(A):
	L = np.linalg.cholesky(A)
	LInv = np.linalg.inv(L)
	AInv = np.dot(LInv.T, LInv)
	logDeterminant = -2.0 * np.sum(np.log(np.diag(LInv)))   # Why the minus sign?
	return AInv, logDeterminant
