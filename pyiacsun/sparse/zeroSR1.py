from __future__ import print_function
__all__ = ["zeroSR1"]

import numpy as np
import scipy.linalg
import datetime

def zeroSR1(fcnGrad, h, prox, options):
	"""
	ZEROSR1 Solves smooth + nonsmooth/constrained optimization problems
	xk,nit, stepSizes = zeroSR1(fcnGrad, h, prox_h, opts)

	This uses the zero-memory SR1 method (quasi-Newton) to solve:

	min_x f(x)  + h(x)

	where
	'fcnGrad' calculates f(x) and its gradient at x,
	and h(x) is a non-smooth term that can be infinite-valued (a constraint),
	so long as you present a function 'prox' that computes diagional plus
	rank-1 projections. The 'prox' function should accept at least three inputs:

	'h' is the non-smooth function, and prox_h is a function with
	3 or 4 inputs that returns:
		y = prox_h( x0 , d, v, ) 
	where
		y = argmin_x h(x) + 1/2||x-x0||^2_B 
	and
		B = inv(H) = inv( diag(D) + v*v' )
	or, for the case with 4 arguments, y = prox_h( x0, d, v, sigma )
	then B = inv( diag(D) + sigma*v*v' ) where sigma should be +1 or -1
	The 4 argument case only matters when opts.SR1=true and opts.BB_type=1
	or opts.SR1=true, opts.BB_type=1 and opts.SR1_diagWeight > 1

	If 'prox_h' isn't provided or is [], it defaults to the identity mapping, which corresponds
	to the case when h=0.

	'prox_h' is mean to be given by something like prox_rank1_l1
	e.g., 
	prox        = @(x0,d,v) prox_rank1_l1( x0, d, v, lambda );
	or, for 4 arguments,
	prox        = @(x0,d,v,varargin) prox_rank1_l1( x0, d, v, lambda, [], varargin{:} );

	"opts" is a dictionary with additional options
	
	opts = {'tol': 1e-6, 'grad_tol' : 1e-6, 'nmax' : 1000, 'N' : N, 'L': normQ, 'verbose': 25}
	
	- 'tol': final tolerance in function
	- 'grad_tol': final tolerance in gradient
	- 'nmax': maximum number of iterations
	- 'N': size of signal (optional)
	- 'x0': initial estimation of the signal (optional)
	- 'L': estimation of the Lipschitz constant (or diagonal scaling)
	- 'verbose': step size for the printing (=0 no printing)
		
	Stephen Becker and Jalal Fadili, Nov 24 2011 -- Dec 2012
	Copied from zeroSR1.m Dec 11 2012
	Feb 28 2014, unnesting all functions to make compatible with octave.

	See also proximalGradient.m	
	
	Python version directly translated from Matlab version (including comments): A. Asensio Ramos (March 12, 2015)
	"""	
	
	start = datetime.datetime.now()
	
	if (('N' in options) & ('x0' not in options)):
		N = options['N']
		xk = np.zeros((N,1))
	elif (('N' not in options) & ('x0' in options)):
		xk = x0.copy()
		N = len(xk)
	else:
		print("I have no way to know the size of the signal to retrieve. Please, set options['N'] or options['x0']")
		sys.exit(1)		
	
	maxStag = 10
	SR1 = True
	BB = True
	nMax = options['nmax']
	
	L = options['L']
	
	Sigma = 1
	BB_type = 2
	if ((SR1) & (BB_type == 1)):
		print("zeroSR1:experimental - With zero-memory SR1, BB_type=1 is an untested feature")
		Sigma = -1
	
	SR1_diagWeight = 0.8*(BB_type==2) + 1.0*(BB_type==1)
	if ((SR1) & (BB_type == 2) & (SR1_diagWeight > 1)):
		Sigma = -1
	
	skipBB = False
	stag = 0
	
	fxOld = np.inf
	t = 1.0 / L
	stepSizes = np.zeros((nMax,1+SR1))
	
# Initialization
	xk_old = xk
	f, gradient = fcnGrad(xk)
	f_xk = np.empty([])
	gradientOld = gradient.copy()
		
# Begin algorithm
	for nIteration in range(nMax):

# "sk" and "yk" are the vectors that will give us quasi-Newton
# information (and also used in BB step, since that can be
# seen as a quasi-Newton method)		
		sk = xk - xk_old
		yk = gradient - gradientOld
		
		if ((nIteration > 0) & (np.linalg.norm(yk) < 1e-13)):
			print("zeroSR1:zeroChangeInGradient. Gradient isn't changing, try changing L")
			yk = np.asarray([])
			skipBB = True
		
# Find and initial stepsize
		if ((BB) & (nIteration > 0) & (not skipBB)):
			if (BB_type == 1):
				t   = np.linalg.norm(sk)**2 / (sk.T.dot(yk)) # eq (1.6) in Dai/Fletcher. This is longer
			else:
				t   = sk.T.dot(yk) / np.linalg.norm(yk)**2 # eq (1.7) in Dai/Fletcher. This is shorter
			if (t < 1e-14):
				print("Curvature condition violated!")
				stag = np.inf			
		
			if (SR1):
# we canot take a full BB step, otherwise we exactly satisfy the secant 
# equation, and there is no need for a rank-1 correction.
				t    = SR1_diagWeight*t # SR1_diagWeights is a scalar less than 1 like 0.6
			
			H0 = lambda x: t*x
			diagH = t*np.ones((N,1))			
		else:
			t = 1.0 / L
			H0 = lambda x: t*x
			diagH = t*np.ones((N,1))
		
		skipBB = False
		stepSizes[nIteration,0] = t
		
				
# ---------------------------------------------------------------------
# -- Quasi-Newton -- Requires: H0, and builds H
# ---------------------------------------------------------------------		
		if ((SR1) & (nIteration > 0) & (yk.size != 0)):
			gs = yk.T.dot(sk)
			if (gs < 0):
				print("Serious curvature condition problem!")
				stag = np.inf
			
			H0 = lambda x: diagH * x
			vk = sk - H0(yk)
			vkyk = vk.T.dot(yk)
			Sigma_local = np.sign(vkyk[0])
			
			if ((Sigma_local * vkyk) <= 0):
				print("Warning: violated curvature conditions")
				vk = []
				H = H0
				stepSizes[nIteration,1] = 0
			else:
				vk /= np.sqrt(Sigma_local * vkyk)
				H = lambda x: H0(x) + Sigma_local * vk.dot(vk.T.dot(x))
				stepSizes[nIteration,1] = vk.T.dot(vk)
		else:
			Sigma_local = Sigma
			H = H0
			vk = []
		
# ---------------------------------
# Make the proximal update
# ---------------------------------
		p = H(-gradient) # Scaled descent direction. H includes the stepsize
		xk_old = xk.copy()
		
		if (Sigma_local != 1):
			xk = prox(xk_old + p, diagH, vk, Sigma_local)
		else:
			xk = prox(xk_old + p, diagH, vk)
			
		norm_grad = np.linalg.norm(xk - xk_old)
		if ( (np.any(np.isnan(xk))) | (np.linalg.norm(xk) > 1.e10)):
			stag = np.inf
			xk = xk_old
			print("Prox algorithm failed, probably due to numerical cancellations")
		
# Update function and gradient
		gradientOld = gradient.copy()
		f_xk, gradient = fcnGrad(xk)
		fx = f_xk + h(xk)
		df = np.abs(fx - fxOld) / np.abs(fxOld)
		fxOld = fx.copy()
		
# Print iteration and test for stopping
		if ((df < options['tol']) | (t < 1e-10) | (np.any(np.isnan(fx))) | (norm_grad < options['grad_tol'])):
			stag += 1
		
		if ((options['verbose'] != 0)):
			if (((nIteration+1) % options['verbose'] == 0) | (stag > maxStag)):
				try:
					print("Iter: {0:5d}, f: {1:.3e}, df: {2:.2e}, ||grad||: {3:.2e}, step: {4:.2e}".format(nIteration+1, fx, df, norm_grad, t[0,0]))
				except:
					print("Iter: {0:5d}".format(nIteration+1))
			
		if (stag > maxStag):
			delta = datetime.datetime.now() - start
			print("Quitting. Reached tolerance. Ellapsed time: {0:2f} s".format(delta.total_seconds()))
			break
		
	return xk, nIteration, stepSizes	