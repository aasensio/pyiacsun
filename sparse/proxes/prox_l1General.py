from __future__ import print_function
import numpy as np
from . import prox_l1
from . import prox_l0

def prox_l1General(x, A, At, lambdaPar, y = None, mu=1.0, verbose=False, threshold='soft'):
		"""
		Proximal projection on the condition |W^* T|_1
		It solves the problem by returning
		sol = argmin_{z} 0.5*||x - z||_2^2 + gamma * ||A z - TRef ||_1
		using the algorithm published in
		M. Fadili and J. Starck. Monotone operator splitting for optimization
		problems in sparse recovery. In Image Processing (ICIP), 2009 16th IEEE
		International Conference on, pages 1461-1464. IEEE, 2009
		
		Args:
		    x (float): input vector
		    A (float): forward operator
		    At (function): backward operator
		    lambdaPar (float): regularization parameter
		    y (function): reference vector that is substracted to A*sol
		    mu (float): norm of the operator
		"""

		u_l1 = np.zeros(A(x).shape)

		if (y == None):
			y = np.zeros(u_l1.shape)
						
		muInv = 1.0 / mu
		sol = x - At(u_l1)
		prev_obj = 0.0
		iteration = 0
		while True:		
			dummy = A(sol)
			
			norm_l1 = lambdaPar * np.sum(np.abs(dummy-y))
			norm_l2 = 0.5 * np.linalg.norm(x - sol, 2)**2 
			norm_obj = norm_l2 + norm_l1			
			if (norm_obj == 0):
				rel_obj = 0.0
			else:
				rel_obj = np.abs(norm_obj - prev_obj) / norm_obj
			
			if (verbose):
				print("  Inner iter {0} - ||Ax||_1={1}, ||x-z||^2={2} - rel_obj={3}".format(iteration, norm_l1, norm_l2, rel_obj))
			
			res = muInv * u_l1 + dummy

			if (threshold == 'soft'):
				dummy = prox_l1(res - y, lambdaPar * muInv) + y
			if (threshold == 'hard'):
				dummy = prox_l0(res - y, lambdaPar * muInv) + y				
			
			u_l1 = 1.0 / muInv * (res - dummy)
			sol = x - At(u_l1)
			
			if (rel_obj < 1e-5):
				break
			if (iteration > 500):
				break
			
			prev_obj = norm_obj
			iteration += 1
		
		return sol