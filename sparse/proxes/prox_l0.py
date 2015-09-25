import numpy as np

def prox_l0(x, lambdaPar):
		"""
		Hard thresholding operator
		"""
		xPar = np.copy(x)		
		xPar[np.abs(x) < lambdaPar] = 0.0
		return xPar