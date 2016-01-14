from __future__ import print_function
import numpy as np
from . import prox_l1
from . import prox_l0

def prox_l2General(x, A, At, lambdaPar, y = None, mu=1.0, verbose=False, threshold='soft', rtol=1e-5, maxIteration=500):
		"""
		Proximal projection on the condition |W^* T|_2
		It solves the problem by returning
		sol = argmin_{z} 0.5*||x - z||_2^2 + gamma * ||A z - TRef ||_2^2
		using FISTA
		
		Args:
		    x (float): input vector
		    A (float): forward operator
		    At (function): backward operator
		    lambdaPar (float): regularization parameter
		    y (function): reference vector that is substracted to A*sol
		    mu (float): norm of the operator
		    verbose (bool, optional): verbose or not
		    threshold (str, optional): 'soft', 'hard' or a function that performs the thresholding. The function
		    							accepts the vector and the thresholding parameter
		"""

		u_n = np.zeros(x.shape)
		sol = np.zeros(x.shape)

		muInv = 1.0 / mu
		
		if (y == None):
			y = np.zeros(x.shape)

		tn = 1.0
		prev_l2 = 0.0
		loop = 0

		stepsize = 1.0 / ((2.0 * lambdaPar)**2 * muInv + 1)
		grad = lambda z : z - x + lambdaPar*2*At(A(z) - y)

		while True:		
			dummy = A(sol - y)
			
			norm_l2 = 0.5 * np.linalg.norm(x - sol, 2)**2 + lambdaPar * np.linalg.norm(dummy)
			if (norm_l2 == 0):
				rel_l2 = 0.0
			else:
				rel_l2 = np.abs(norm_l2 - prev_l2) / norm_l2
			
			if (verbose):
				print("  Inner iter {0} - ||Ax-y||_2^2={1} - rel_obj={2}".format(loop, norm_l2, rel_l2))
			
			if (rel_l2 < rtol):
				break
			if (loop > maxIteration):
				break
			
			x_n = u_n - stepsize * grad(u_n)
			tn1 = 0.5 * (1.0 + np.sqrt(1.0 + 4.0 * tn**2))
			u_n = x_n + (tn-1) / tn1 * (x_n - sol)

			sol = np.copy(x_n)
			tn = np.copy(tn1)

			prev_l2 = np.copy(norm_l2)
			loop += 1

		return sol