import numpy as np

def prox_xlim(x, left, right, verbose=False):
		"""
		Proximal projection on a box
		It solves the problem by returning
		sol = argmin_{z} 0.5*||x - z||_2^2, s.t. x < right and x > left
		
		Args:
		    x (float): input vector
		    left (TYPE): left bound
		    right (TYPE): right bound
		    verbose (bool, optional): verbose
				
		"""
		sol = np.copy(x)
		sol[x > right] = right
		sol[x < left] = left

		return sol