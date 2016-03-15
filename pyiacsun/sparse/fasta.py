from __future__ import print_function

__all__ = ["fasta"]

import numpy as np
import scipy.linalg
import datetime

class fasta(object):

	def __init__(self, A, At, f, gradf, g, proxg, x0, maxIter=1000, tol=1e-3, verbose=False, recordObjective=False,
		recordIterates=False, adaptive=True, accelerate=False, restart=True, backtrack=True, stepsizeShrink=0.2,
		window=10, eps_r=1e-8, eps_n=1e-8, L=0.0, tau=0.0):
		"""
		This method solves the problem
	        minimize f(Ax)+g(x)
	   Where A is a matrix, f is differentiable, and both f and g are convex.  
	   The algorithm is an adaptive/accelerated forward-backward splitting.  
	   The user supplies function handles that evaluate 'f' and 'g'.  The user 
	   also supplies a function that evaluates the gradient of 'f' and the
	   proximal operator of 'g', which is given by
	                proxg(z,t) = argmin t*g(x)+.5||x-z||^2.
	
	  Inputs:
	    A     : 
	    At    : The adjoint (transpose) of 'A.' Optionally, a function handle
	             may be passed.  
	    gradf : A function of z, computes the gradient of f at z
	    proxg : A function of z and t, the proximal operator of g with
	             stepsize t.  
			
	Args:
	    A (TYPE): a function that returns A*x
	    At (TYPE): a function that returns A'*x
	    f (TYPE): a function of x that computes f
	    gradf (TYPE): a function of x that computes the gradient of f at x
	    g (TYPE): a function of x that computes g
	    proxg (TYPE): a function of x and t that computes the proximal operator of g with step size t
	    x0 (TYPE): initial guess	   
	"""

		self.A = A
		self.At = At
		self.f = f
		self.gradf = gradf
		self.g = g
		self.proxg = proxg
		self.x0 = x0

# Options
		self.maxIter = maxIter
		self.tol = tol
		self.verbose = verbose
		self.recordObjective = recordObjective
		self.recordIterates = recordIterates
		self.adaptive = adaptive
		self.accelerate = accelerate
		self.restart = restart
		self.backtrack = backtrack
		self.stepsizeShrink = stepsizeShrink
		self.window = window
		self.eps_r = eps_r
		self.eps_n = eps_n
		self.L = L
		self.tau = tau

	# The adaptive method can expand the stepsize, so we choose an aggressive value by default
	# If the stepsize is monotonically decreasing, we don't want to make it smaller than we need
		if (not self.adaptive or self.accelerate):
			self.stepsizeShrink = 0.5

	# Define the algorithm used
		self.mode = 'plain'
		if (self.adaptive):
			self.mode = 'adaptive'
		if (self.accelerate):
			if (self.restart):
				self.mode = 'accelerated(FISTA)+restart'
			else:
				self.mode = 'accelerated(FISTA)'

		if (self.L <= 0.0 and self.tau <= 0.0):
			x1 = np.random.normal(size=self.x0.shape)
			x2 = np.random.normal(size=self.x0.shape)
			gradf1 = self.At(self.gradf(self.A(x1)))
			gradf2 = self.At(self.gradf(self.A(x2)))
			self.L = np.linalg.norm(gradf1 - gradf2) / np.linalg.norm(x2 - x1)
			self.L = np.max([self.L, 1e-6])
			self.tau = 2.0 / self.L / 10

		if (self.tau == 0 and self.L != 0):
			self.tau = 1.0 / self.L

		if (self.L == 0 and self.tau != 0):
			self.L = 1.0 / self.tau

	def checkAdjoint(self):
		"""Check if A and At are adjoint operators
			
		Returns:
		    TYPE: error in the two inner products
		"""
		xTest = np.random.normal(size=self.x0.shape)
		Ax = self.A(xTest)
		yTest = np.random.normal(size=Ax.shape)
		Aty = self.At(yTest)

		inner1 = Ax.dot(yTest)
		inner2 = xTest.T.dot(Aty)

		error = np.abs(inner1 - inner2) / np.max([np.abs(inner1),np.abs(inner2)])

		return error

	def optimize(self):
		"""
			This method solves the problem
		        minimize f(Ax)+g(x)
		   Where A is a matrix, f is differentiable, and both f and g are convex.  
		   The algorithm is an adaptive/accelerated forward-backward splitting.  
		   The user supplies function handles that evaluate 'f' and 'g'.  The user 
		   also supplies a function that evaluates the gradient of 'f' and the
		   proximal operator of 'g', which is given by
		                proxg(z,t) = argmin t*g(x)+.5||x-z||^2.
			  	
		Args:
		    A (TYPE): a function that returns A*x
		    At (TYPE): a function that returns A'*x
		    f (TYPE): a function of x that computes f
		    gradf (TYPE): a function of x that computes the gradient of f at x
		    g (TYPE): a function of x that computes g
		    proxg (TYPE): a function of x and t that computes the proximal operator of g with step size t
		    x0 (TYPE): initial guess
		    options (TYPE): dictionary
		"""
		
		x0 = np.copy(self.x0)
		tau1 = np.copy(self.tau)
		error = self.checkAdjoint()
		if (error > 1e-9):
			print ("At is not the adjoint of A")

		if (self.verbose):
			print("FASTA: {0} - Max. Iter: {1}".format(self.mode, self.maxIter))

		residual = np.zeros(self.maxIter)
		normalizedResidual = np.zeros(self.maxIter)
		taus = np.zeros(self.maxIter)
		fVals = np.zeros(self.maxIter)               # The value of 'f', the smooth objective term
		objective = np.zeros(self.maxIter)           # The value of the objective function (f+g)
		funcValues = np.zeros(self.maxIter)          # Values of the optional 'function' argument in 'opts'
		totalBacktracks = 0                          # How many times was backtracking activated?
		backtrackCount = 0                           # Backtracks on this iterations

# Initialize array values
		x1 = np.copy(x0)
		d1 = self.A(x1)
		f1 = self.f(d1)
		fVals[0] = f1
		gradf1 = self.At(self.gradf(d1))

		if (self.accelerate):
			xAccel1 = np.copy(x0)
			dAccel1 = np.copy(d1)
			alpha1 = 1.0

		maxResidual = -np.inf
		minObjectiveValue = np.inf

		if (self.recordObjective):
			objective[0] = f1 + self.g(x0)

		start = datetime.datetime.now()

		for i in range(self.maxIter):

# Rename iterates relative to loop index.  "0" denotes index i, and "1" denotes index i+1
			x0 = np.copy(x1)
			gradf0 = np.copy(gradf1)
			tau0 = np.copy(tau1)

# FBS step: obtain x_{i+1} from x_i
			x1hat = x0 - tau0*gradf0
			x1 = self.proxg(x1hat, tau0)

# Non-monotone backtracking line search
			Dx = x1 - x0
			d1 = self.A(x1)
			f1 = self.f(d1)
			if (self.backtrack):
				left = np.max([i-self.window, 0])
				right = np.max([i-1, 0])

				if (left == right):
					M = fVals[0]
				else:
					M = np.max(fVals[left:right])

				backtrackCount = 0

				temp = Dx.dot(gradf0) + np.linalg.norm(Dx, 2)**2 / (2.0*tau0)

				while ( ( (f1 - 1e-12) > (M + temp) ) and (backtrackCount < 10)):
					tau0 *= self.stepsizeShrink
					x1hat = x0 - tau0 * gradf0
					x1 = self.proxg(x1hat, tau0)
					d1 = self.A(x1)
					f1 = self.f(d1)
					Dx = x1 - x0
					backtrackCount += 1
					temp = Dx.dot(gradf0) + np.linalg.norm(Dx, 2)**2 / (2.0*tau0)					

				totalBacktracks += backtrackCount

			if (self.verbose and backtrackCount > 10):
				print ("WARNING: excessive backtracking ({0} steps) - current step is {1}".format(backtrackCount, tau0))

# Record convergence information
			taus[i] = tau0
			residual[i] = np.linalg.norm(Dx) / tau0
			maxResidual = np.max([maxResidual, residual[i]])
			normalizer = np.max([np.linalg.norm(gradf0), np.linalg.norm((x1 - x1hat)/tau0)]) + self.eps_n
			normalizedResidual[i] = residual[i] / normalizer
			fVals[i] = f1

			newObjectiveValue = residual[i]

# Method is non-monotonic. Record the best solution
			if (newObjectiveValue < minObjectiveValue):
				bestObjectiveIterate = x1
				minObjectiveValue = newObjectiveValue

			if (self.verbose):
				if (i % 100 == 0):					
					print("It: {0:6d} - resid: {1:12.4e} - backtrack: {2:2d} - tau: {3:12.4e}".format(i, residual[i], backtrackCount, float(tau0)))
			
# Check convergence
			if ((residual[i] / (maxResidual + self.eps_r) < self.tol) or (normalizedResidual[i] < self.tol)):				
				return bestObjectiveIterate

# Adaptive
			if (self.adaptive and not self.accelerate):
				gradf1 = self.At(self.gradf(d1))
				Dg = gradf1 + (x1hat - x0) / tau0
				dotprod = Dx.dot(Dg)
				tau_s = np.linalg.norm(Dx)**2 / dotprod
				tau_m = dotprod / np.linalg.norm(Dg)**2
				tau_m = np.max([tau_m, 0])
				if (2*tau_m > tau_s):
					tau1 = tau_m
				else:
					tau1 = tau_s - 0.5 * tau_m
				if (tau1 <= 0 or np.isinf(tau1) or np.isnan(tau1)):
					tau1 = 1.5 * tau0

# Acceleration
			if (self.accelerate):
				xAccel0 = np.copy(xAccel1)
				dAccel0 = np.copy(dAccel1)
				alpha0 = np.copy(alpha1)
				xAccel1 = np.copy(x1)
				dAccel1 = np.copy(d1)

				# Check for restarting
				if (self.restart and np.dot(x0-x1, x1-xAccel0) > 0):
					alpha0 = 1.0

				alpha1 = 0.5 * (1.0 + np.sqrt(1.0 + 4.0 * alpha0**2))

				x1 = xAccel1 + (alpha0-1.0) / alpha1 * (xAccel1 - xAccel0)
				d1 = dAccel1 + (alpha0-1.0) / alpha1 * (dAccel1 - dAccel0)

				gradf1 = self.At(self.gradf(d1))
				fVals[i] = self.f(d1)
				tau1 = tau0

			if (not self.adaptive and not self.accelerate):
				gradf1 = self.At(self.gradf(d1))
				tau1 = np.copy(tau0)

		return bestObjectiveIterate
