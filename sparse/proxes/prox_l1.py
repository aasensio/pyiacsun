__all__ = ["prox_l1"]

import numpy as np

def prox_l1(x, lambdaPar):
	"""
	Soft thresholding operator
	
	Args:
	    x (TYPE): Description
	    lambdaPar (TYPE): Description
	"""
	return np.sign(x) * np.fmax(np.abs(x) - lambdaPar, 0)