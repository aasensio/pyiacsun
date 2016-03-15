__all__ = ["normSquared"]

import numpy as np

def normSquared(x, A, At, b=0, c=0, constant=0):
	"""
	f, g = normSquared(x,A,At,b,c,constant)
	returns the objective function f and its gradient
	
	A is a function that computes the matrix-vector product
	At is a function that computes the matrix-vector product for the transposed matrix
	
	By default, b=0 and c=0
		
	March 4 2014, Stephen Becker, stephen.beckr@gmail.com
	
	See also quadraticFunction.m
	
	Python version: A. Asensio Ramos (March 12, 2015)
	
	"""
	res = A(x) - b
	cx = 0.0
	if (c != 0):
		cx = np.dot(c, x)
		
	f = 0.5 * np.linalg.norm(res)**2 + cx + constant
	g   = At(res) + c
	
	return f, g
