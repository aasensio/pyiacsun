from __future__ import print_function

__all__ = ['rvm']
import numpy as np

class rvm(object):
	
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Constructor
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	def __init__(self, basis, targets, noise=None):		
		"""
		Class that implements a Relevance Vector Machine (RVM) as described in Tipping, M. E. and A. C. Faul (2003)
		The RVM does regression using a linear function with a potentially very large N functions using
		a sparsity constraint on the vector w_i
		
		y(x) = w_1*phi_1(x)+w_2*phi_2(x)+...+w_N*phi_N(x)
		
		Instantiate the class using
		
		p = rvm(Basis, Outputs, noise=0.018)
		
		 - Outputs: value of the function at nPoints locations x_i
		 - Basis: matrix of size [nPoints x nFunctions] of the nFunctions evaluated at the x_i locations of the observations
		 - noise: optional parameter to give an estimation of the noise standard deviation. It can be absent and the noise variance
		          will be estimated
		          
		The RVM is fitted using
		
		p.iterateUntilConvergence()
		
		or individual iterations can be carried out using 
		
		p.iteration()
		
		The value of the weights w_i is returned in p.wInferred, so that the function evaluated at the points is
		given by
		
		np.dot(Basis, p.wInferred)
		
		You can also evaluate the Basis functions in a finer x_i grid to produce a smoother regressor
		
		"""
		self.basis = 1.0*basis
		self.targets = targets
		self.noise = noise

# A "reasonable" initial value for the noise in the Gaussian case
		self.gaussianSNRInit = 0.1
		
# "Reasonable" initial alpha bounds
		self.initAlphaMax = 1.0e3
		self.initAlphaMin = 1.0e-3
		
# If noise standard deviation is given, compute beta from the noise standard deviation
		if (noise is not None):
			self.beta = 1.0 / noise**2
			
# Initialize basis (phi), mu and alpha

# First compute linearised output for use in heuristic initialization
		self.TargetsPseudoLinear = self.targets
		
# BASIS PREPROCESSING:
# Scale basis vectors to unit norm. This eases some calculations and
# will improve numerical robustness later.
		self.preprocessBasis()
				
		self.initialization()
		
# Cache some quantities
		self.basisPHI = np.dot(self.basis.T, self.PHI)
		self.basisTargets = np.dot(self.basis.T, self.targets)

		
# Full computation
#
# Initialise with a full explicit computation of the statistics
#
# NOTE: The AISTATS paper uses "S/Q" (upper case) to denote the key
# "sparsity/quality" Factors for "included" basis functions, and "s/q"
# (lower case) for the factors calculated when the relevant basis
# functions are "excluded".
#
# Here, for greater clarity:
#
#	All S/Q are denoted by vectors S_in, Q_in
#	All s/q are denoted by vectors S_out, Q_out

		self.fullStatistics()
		
		self.N, self.MFull = self.basis.shape
		self.M = self.PHI.shape[1]
		
		self.addCount = 0
		self.deleteCount = 0
		self.updateCount = 0
		self.maxLogLike = -1.e10
		
		self.controls = {'BetaUpdateStart': 10, 'BetaUpdateFrequency': 5, 'ZeroFactor': 1.e-12, 'PriorityAddition': 0, 
			'PriorityDeletion': 1, 'AlignmentMax': 0.999e0, 'MinDeltaLogAlpha': 1.e-3}
		
		self.logMarginalLog = np.asarray([self.logML])
		self.countLoop = 0
		self.alignDeferCount = 0
		self.loop = 0
		self.lastIteration = False
		self.actionReestimate = 0
		self.actionAdd = 1
		self.actionDelete = -1
		self.actionTerminate = 10
		self.actionNoiseOnly = 11
		self.actionAlignmentSkip = 12
		self.selectedAction = 0
		self.alignedOut = []
		self.alignedIn = []

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Initialization
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	def initialization(self):
		
# 1) the starting basis, PHI

# Take account of "free basis": it needs to be included from the outset
# Set initial basis to be the largest projection with the targets		
		proj = np.dot(self.basis.T, self.TargetsPseudoLinear)
		
		self.Used = np.atleast_1d(np.argmax(abs(proj)))
		foo = proj[self.Used]
		
		self.PHI = np.atleast_2d(self.basis[:,self.Used])
		M = self.Used.shape
		
# 2) the most probable weights
#  mu will be calculated analytically later in the Gaussian case
# 	mu = []

# 3) the hyperparameters

# Exact for single basis function case (diag irrelevant),
# heuristic in the multiple case	
		p = np.diag(np.dot(self.PHI.T, self.PHI)) * self.beta
		q = np.dot(self.PHI.T, self.targets) * self.beta
		self.Alpha = np.asarray([p**2 /(q**2-p)])
	
# Its possible that one of the basis functions could already be irrelevant (alpha<0), so trap that		
		ind = np.where(self.Alpha < 0.0)[0]
		if (ind.shape[0] != 0):
			self.Alpha[ind] = self.initAlphaMax
	
		print("Initial alpha = {0}".format(self.Alpha))

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Pre-process basis
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	def preprocessBasis(self):
		
		N, M = self.basis.shape
		
		self.scales = np.sqrt(np.sum(self.basis**2,0))
		self.scales[self.scales == 0] = 0.0
						
		for i in range(M):
			self.basis[:,i] = self.basis[:,i] / self.scales[i]
			
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Compute the full statistics
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	def fullStatistics(self):
		N, MFull = self.basis.shape
		M = self.PHI.shape[1]
		
# Cholesky decomposition to carry out the inverse and compute Sigma matrix
		
		U = np.linalg.cholesky(np.dot(self.PHI.T,self.PHI) * self.beta + np.diag(self.Alpha))
		Ui = np.linalg.inv(U)
		#pdb.set_trace()
		self.Sigma = np.dot(Ui, Ui.T)
		
# Posterior mean
		self.Mu = np.dot(self.Sigma, np.dot(self.PHI.T, self.targets)) * self.beta
		
# Data error and likelihood		
		y = np.dot(self.PHI, self.Mu)
		e = (self.targets - y)
		ED = np.dot(e.T, e)
		dataLikelihood = 0.5 * (N * np.log(self.beta) - self.beta * ED)		
		
# Compute the log-marginal posterior
		logDetHOver2 = np.sum(np.log(np.diag(U)))
		self.logML = dataLikelihood - 0.5 * np.dot(self.Mu.T**2, self.Alpha) + 0.5 * np.sum(np.log(self.Alpha)) - logDetHOver2
		
# Well-determinedness factors
		diagC = np.sum(Ui**2, 1)
		self.Gamm = 1.0 - self.Alpha * diagC
		#pdb.set_trace()
		
# Compute Q & S values
# Q: "quality" factor - related to how well the basis function contributes
# to reducing the error
# S: "sparsity factor" - related to how orthogonal a given basis function
# is to the currently used set of basis functions
		self.betaBasisPHI = np.dot(self.basis.T, self.PHI * self.beta * np.ones((1,M)))
		self.SIn = self.beta - np.sum(np.dot(self.betaBasisPHI,Ui)**2,1,keepdims=True)
		self.QIn = np.atleast_2d(self.beta * (self.basisTargets - np.dot(self.basisPHI, self.Mu))).T
		
		self.SOut = 1*self.SIn
		self.QOut = 1*self.QIn
		
		self.SOut[self.Used] = (self.Alpha * self.SIn[self.Used]) / (self.Alpha - self.SIn[self.Used])
		self.QOut[self.Used] = (self.Alpha * self.QIn[self.Used]) / (self.Alpha - self.SIn[self.Used])
				
		self.factor = (self.QOut**2 - self.SOut)
		
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# One iteration
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	def oneIteration(self):

#*****************
# Decision phase
#*****************

# Compute change in likelihood
		deltaML = np.zeros((self.MFull, 1))
		action = self.actionReestimate * np.ones((self.MFull,1))
		usedFactor = self.factor[self.Used]
		N, M = self.PHI.shape
		
		
# Reestimation
		iu = np.where(usedFactor > self.controls['ZeroFactor'])[0]
		index = self.Used[iu]
		newAlpha = self.SOut[index]**2 / self.factor[index]
		delta = 1.0 / newAlpha - 1.0 / self.Alpha[iu]
		
# Quick computation of change in log-likelihood given all re-estimations
		deltaML[index] = 0.5 * (delta * self.QIn[index]**2 / (delta * self.SIn[index]+1.0) - np.log(1.0+self.SIn[index]*delta))
		
# Deletion
		iu = np.where(usedFactor <= self.controls['ZeroFactor'])[0]
		index = self.Used[iu]
		anyToDelete = (not index.size == 0) and M > 1
		
		if (anyToDelete):
			deltaML[index] = -0.5 * (self.QOut[index]**2 / (self.SOut[index] + self.Alpha[iu]) - np.log(1.0+self.SOut[index] / self.Alpha[iu]))
			action[index] = self.actionDelete
			
# Addition
		GoodFactor = self.factor > self.controls['ZeroFactor']
		GoodFactor[self.Used] = 0
		GoodFactor[self.alignedOut] = 0
		index = np.where(GoodFactor)[0]
		anyToAdd = (not index.size == 0)
		
		if (anyToAdd):
# Quick computation of change in log-likelihood
			quot = self.QIn[index]**2 / self.SIn[index]
			deltaML[index] = 0.5 * (quot - 1.0 - np.log(quot))
			action[index] = self.actionAdd
			
		if (anyToAdd and self.controls['PriorityAddition']) or (anyToDelete and self.controls['PriorityDeletion']):
			
# We won't perform re-estimation this iteration
			deltaML[action == self.actionReestimate] = 0
			
# We should enforce add if preferred and delete
			if (anyToAdd and self.controls['PriorityAddition'] and (not self.controls['PriorityDeletion'])):
				deltaML[action == self.actionDelete] = 0
				
			if (anyToDelete and self.controls['PriorityDeletion'] and (not self.controls['PriorityAddition'])):
				deltaML[action == self.actionAdd] = 0

# Finally, choose the action that results in the greatest change in likelihood
		nu = np.atleast_1d(np.argmax(deltaML))	
		self.deltaLogMarginal = deltaML[nu]
		self.selectedAction = action[nu]
		anyWorthWhileAction = self.deltaLogMarginal > 0
		
# If basis nu is already in the model, find its index
		j = []
		if (self.selectedAction == self.actionReestimate) or (self.selectedAction == self.actionDelete):
			j = np.where(self.Used == nu)[0]
			
		self.Phi = np.atleast_2d(self.basis[:,nu])
		newAlpha = self.SOut[nu]**2 / self.factor[nu]
		
		change = np.abs(np.log(newAlpha) - np.log(self.Alpha[j]))
							
		if (not anyWorthWhileAction) or ((self.selectedAction == self.actionReestimate) and (change < self.controls['MinDeltaLogAlpha']) and (not anyToDelete)):
			self.selectedAction = self.actionTerminate
			
# Alignment checks for addition		
		if (self.selectedAction == self.actionAdd):
			p = np.dot(self.Phi.T, self.PHI)
			findAligned = np.where(p > self.controls['AlignmentMax'])[0]
			numAligned = findAligned.size
			
			if (numAligned > 0):
# The added basis function is effectively indistinguishable from one present already
				self.selectedAction = self.actionAlignmentSkip
				self.alignDeferCount += 1
				
# Take note not to try this again
				self.alignedOut = np.append(self.alignedOut, nu*np.ones((numAligned,1))).astype(int)
				self.alignedIn = np.append(self.alignedIn, self.Used[findAligned])
		
# Alignment checks for deletion
		if (self.selectedAction == self.actionDelete):
			findAligned = np.where(self.alignedIn == nu)[0]
			numAligned = findAligned.size
			if (numAligned > 0):
				reinstated = self.alignedOut[findAligned]
				self.alignedIn = np.delete(self.alignedIn, findAligned, 0)
				self.alignedOut = np.delete(self.alignedOut, findAligned, 0)
				
#*****************
# Action phase
#*****************
		updateRequired = False
		
		if (self.selectedAction == self.actionReestimate):			
# Basis function nu is already in the model and we're reeestimatig its alpha
			oldAlpha = self.Alpha[j]
			self.Alpha[j] = newAlpha
			S_j = self.Sigma[:,j]
			deltaInv = 1.0 / (newAlpha - oldAlpha)
			kappa = 1.0 / (self.Sigma[j,j] + deltaInv)
			tmp = kappa * S_j
			newSigma = self.Sigma - np.dot(tmp, S_j.T)
			deltaMu = -self.Mu[j] * tmp
			self.Mu += deltaMu
			
			self.SIn += kappa * np.dot(self.betaBasisPHI, S_j)**2
			self.QIn -= np.dot(self.betaBasisPHI, deltaMu)
			
			self.updateCount += 1
			updateRequired = True
			
		elif (self.selectedAction == self.actionAdd):						
			self.basisPhi = np.dot(self.basis.T, self.Phi)
			self.basisPHI = np.hstack((self.basisPHI, self.basisPhi))
			self.BPhi = self.beta * self.Phi
			self.BASISBPhi = self.beta * self.basisPhi
			
			tmp = np.dot(np.dot(self.BPhi.T, self.PHI), self.Sigma).T
			
			self.Alpha = np.vstack((self.Alpha,newAlpha))
			self.PHI = np.hstack((self.PHI, self.Phi))
			s_ii = 1.0 / (newAlpha + self.SIn[nu])
			s_i = -s_ii * tmp
			TAU = -np.dot(s_i, tmp.T)
			
			t1 = np.hstack((self.Sigma+TAU, s_i))
			t2 = np.hstack((s_i.T, s_ii))
			newSigma = np.vstack((t1,t2))
			mu_i = s_ii * self.QIn[nu]
			deltaMu = np.vstack((-mu_i*tmp, mu_i))
			self.Mu = np.vstack((self.Mu, 0)) + deltaMu
			
			mCi = self.BASISBPhi - np.dot(self.betaBasisPHI, tmp)
			self.SIn -= s_ii * mCi**2
			self.QIn -= mu_i * mCi
			
			self.Used = np.hstack((self.Used, nu))
			self.addCount += 1
			updateRequired = True						
			
		elif (self.selectedAction == self.actionDelete):			
			self.basisPHI = np.delete(self.basisPHI, j, 1)
			self.PHI = np.delete(self.PHI, j, 1)
			self.Alpha = np.delete(self.Alpha, j, 0)
			s_jj = self.Sigma[j,j]
			s_j = self.Sigma[:,j]
			tmp = s_j / s_jj
			newSigma = self.Sigma - np.dot(tmp, s_j.T)
			newSigma = np.delete(newSigma, j, 0)
			newSigma = np.delete(newSigma, j, 1)
			deltaMu = -self.Mu[j] * tmp
			mu_j = self.Mu[j]
			self.Mu += deltaMu
			self.Mu = np.delete(self.Mu, j, 0)
			
			jPm = np.dot(self.betaBasisPHI, s_j)
			self.SIn += jPm**2 / s_jj
			self.QIn += jPm * mu_j / s_jj
			
			self.Used = np.delete(self.Used, j, 0)
			self.deleteCount += 1
			updateRequired = True
			
		M = len(self.Used)
		
		if (updateRequired):
			self.SOut[:] = self.SIn
			self.QOut[:] = self.QIn
			
			tmp = self.Alpha / (self.Alpha - self.SIn[self.Used])
			
			self.SOut[self.Used] = tmp * self.SIn[self.Used]
			self.QOut[self.Used] = tmp * self.QIn[self.Used]
						
			self.factor = (self.QOut * self.QOut - self.SOut)
			self.Sigma = newSigma
			self.Gamm = 1.0 - self.Alpha * np.atleast_2d(np.diag(self.Sigma)).T
			
			self.betaBasisPHI = self.beta * self.basisPHI
			
			self.logML += self.deltaLogMarginal[0,0]
			self.countLoop += 1
			self.logMarginalLog = np.append(self.logMarginalLog, self.logML)
									
# Something went wrong. Recompute statistics
		if (np.sum(self.Gamm) < 0):
			self.fullStatistics()
		
# Recompute noise if not given
		if (self.noise is None) and ((self.selectedAction == self.actionTerminate) or (self.loop <= self.controls['BetaUpdateStart']) or 
			(self.loop % self.controls['BetaUpdateFrequency'] == 0)):
# Gaussian noise estimate
			betaZ1 = beta
			y = np.dot(self.PHI, self.Mu)
			e = self.targets - y
			beta = (N - np.sum(self.Gamm)) / np.dot(e.T, e)
			
# Work-around zero-noise issue
			beta = np.amin([beta, 1.e6 / np.variance(self.targets)])
			deltaLogBeta = np.log(beta) - np.log(betaZ1)
			
# Full re-computation of statistics after beta update
			if (np.abs(deltaLogBeta) > 1.e-6):
				self.fullStatistics()
				self.countLoop += 1
				self.logMarginalLog = np.append(self.logMarginalLog, self.logML)
				
				if (self.selectedAction == self.actionTerminate):
					self.selectedAction = self.actionNoiseOnly
					print("Noise update. Termination deferred")
		
		#self.AInv = np.diag(1.0/self.Alpha[:,0])
		#Sigma = 1.0/self.beta * np.identity(self.N) + np.dot(np.dot(self.PHI, self.AInv), self.PHI.T)
		#CInv, logD = cholInvert(Sigma)
		#logL = -0.5*logD - 0.5*np.dot(np.dot(self.targets.T,CInv),self.targets)
		
		print("{0:4d} - L={1:10.7f} - Gamma={2:10.7f} (M={3:4d}) - s={4:6.4f}".format(self.loop,self.logML[0]/N, np.sum(self.Gamm), M, np.sqrt(1.0/self.beta)))
				
		
		if (self.selectedAction == self.actionTerminate):
			print("Stopping at iteration {0} - max_delta_ml={1}".format(self.loop, self.deltaLogMarginal[0,0]))
			print("L={0} - Gamma={1} (M={2}) - s={3}".format(self.logML[0]/N, np.sum(self.Gamm), M, np.sqrt(1.0/self.beta)))
		
		iterationLimit = self.loop == 200
		self.lastIteration = iterationLimit
		
											
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Iterate until convergence
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	def iterateUntilConvergence(self):
		while (not self.lastIteration and self.selectedAction != self.actionTerminate):
			self.loop += 1
			self.oneIteration()
		self.postProcess()
			
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Post-process
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	def postProcess(self):		
		#self.fullStatistics()				
		self.index = np.argsort(self.Used)
		self.relevant = self.Used[self.index]
		self.value = self.Mu[self.index] / np.atleast_2d(self.scales[self.Used[self.index]]).T
		self.Alpha = self.Alpha[self.index] / (np.atleast_2d(self.scales[self.Used[self.index]]).T)**2
				
		self.wInferred = np.zeros((self.MFull,1))
		self.wInferred[self.relevant] = self.value
		
#***************
# Test
#***************
#import matplotlib.pyplot as pl
#N = 10
#M = 11
#XStenflo = np.asarray([-2.83000,-1.18000,0.870000,1.90000,3.96000,5.01000,6.25000,8.10000,9.98000,12.1200]).T
#Outputs = np.asarray([0.0211433,0.0167467,0.00938793,0.0183543,-0.00285475,-0.000381000,0.00374350,0.000126900,0.0121750,0.0268133]).T

## Define the basis functions. Gaussians of width 3.4 evaluated at the observed points
#basisWidth = 3.4
#C = XStenflo[:,np.newaxis]
#Basis = np.exp(-(XStenflo-C)**2 / basisWidth**2)
#Basis = np.hstack([Basis, np.ones((1,N)).T])

## Instantitate the RVM object and train it
#p = rvm(Basis, Outputs, noise=0.018)
#p.iterateUntilConvergence()

## Do some plots
#f = pl.figure(num=0)
#ax = f.add_subplot(1,1,1)
#ax.plot(XStenflo, Outputs, 'ro')
#ax.plot(XStenflo, np.dot(Basis, p.wInferred))
