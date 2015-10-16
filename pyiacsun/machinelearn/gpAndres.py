__all__ = ['GaussianProcess']

import numpy as np
import matplotlib.pyplot as pl
from ..linalg import cholInvert
import numpy.core.umath_tests
import scipy.optimize as opt

class GaussianProcess:
	
	def __init__(self, xInput, lambdaGPInitial=1.0, sigmaGPInitial=1.0):
		
		self.lambdaGP = lambdaGPInitial
		self.sigmaGP = sigmaGPInitial
		self.xInput = xInput
		self.CInv = 0
		
	def covariance(self, x, pars):
		"""
		Covariance matrix. It returns the value of the covariance matrix
		and the derivative of the matrix wrt to the hyperparameters
		"""
		lambdaGP = pars[0]
		sigmaGP = pars[1]
		
		expon = np.exp(-0.5 * lambdaGP * x**2)
		
# Covariance matrix
		K = sigmaGP * expon
		
# Derivatives of the covariance matrix
		dKdsigmaGP = expon
		dKdlambdaGP = -0.5 * K * x**2
		
		return K, dKdlambdaGP, dKdsigmaGP

# Returns the marginal likelihood for a Gaussian process
	def marginalLikelihood(self, pars, *args):
		xInput = self.xInput
		yInput, sigmaNoise = args
			
		lambdaGP = np.exp(pars[0])
		sigmaGP = np.exp(pars[1])
								
		K, dKdl, dKds = self.covariance(xInput[np.newaxis,:]-xInput[:,np.newaxis], [lambdaGP, sigmaGP])
				
		C = K + sigmaNoise**2 * np.identity(len(xInput))
		
		CInv, logD = cholInvert(C)
		likelihood = 0.5 * np.dot(np.dot(yInput.T,CInv),yInput) + 0.5 * logD
		
# Jacobian
		jac = np.zeros(2)
		
# dLdlambda
		yInput2 = np.dot(CInv, yInput)
		jac[0] = -0.5 * np.sum(numpy.core.umath_tests.inner1d(CInv, dKdl.T)) + 0.5*np.dot(np.dot(yInput2.T,dKdl),yInput2)
		jac[0] = -jac[0] * lambdaGP
		
# dLdsigma
		jac[1] = -0.5 * np.sum(numpy.core.umath_tests.inner1d(CInv, dKds.T)) + 0.5*np.dot(np.dot(yInput2.T,dKds),yInput2)
		jac[1] = -jac[1] * sigmaGP
					
		return likelihood, jac
	
	def fit(self, y, sigmaNoise):
		x0 = [np.log(self.lambdaGP), np.log(self.sigmaGP)]		
		res = opt.minimize(self.marginalLikelihood, x0, method='BFGS', jac=True, args=(y, sigmaNoise))
		self.lambdaGP, self.sigmaGP = np.exp(res.x)
		self.yNoise = y
		
# Final evaluation of the covariance matrix
		K, dKdl, dKds = self.covariance(self.xInput[np.newaxis,:]-self.xInput[:,np.newaxis], [self.lambdaGP, self.sigmaGP])
				
		C = K + sigmaNoise**2 * np.identity(len(self.xInput))

		self.CInv, logD = cholInvert(C)
		
	def predict(self, xStar):
		K_XStar_XStar, _, _ = self.covariance(xStar[np.newaxis,:]-xStar[:,np.newaxis], [self.lambdaGP, self.sigmaGP])
		K_XStar_X, _, _ = self.covariance(xStar[:,np.newaxis]-self.xInput[np.newaxis,:], [self.lambdaGP, self.sigmaGP])
		K_X_XStar, _, _ = self.covariance(self.xInput[:,np.newaxis]-xStar[np.newaxis,:], [self.lambdaGP, self.sigmaGP])
		
		predicted = K_XStar_X.dot(self.CInv).dot(self.yNoise)
		covariance = K_XStar_XStar - K_XStar_X.dot(self.CInv).dot(K_X_XStar)
		
		return predicted, covariance

if __name__ == "__main__":
	
# Example with a simple case
	nPoints = 50
	sigmaNoise = 0.7
	x = np.linspace(0,8,nPoints)
	y = np.sin(x)
	yNoise = y + sigmaNoise*np.random.randn(nPoints)

	gp = GaussianProcess(x)
	gp.fit(y, sigmaNoise)
	res, variance = gp.predict(x)

	pl.plot(x, yNoise)
	pl.plot(x, res)
	pl.plot(x, y)
