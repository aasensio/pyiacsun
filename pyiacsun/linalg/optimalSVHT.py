from __future__ import print_function
__all__ = ['optimalSVHT']
import numpy as np

def optimalSVHT(matrix):
	"""
	Returns the optimal threshold of the singular values of a matrix, i.e., it computes the optimal number
	of singular values to keep or the optimal rank of the matrix.
	Given a matrix Y=X+sigma*Z, with Y the observed matrix, X the matrix that we want to find, Z an matrix whose elements are iid Gaussian
	random numbers with zero mean and unit variance and sigma is the standard deviation of the noise, the matrix X can be recovered in the
	least squares sense using a truncated SVD (compute the SVD of the matrix, keep only the relevant singular values and reconstruct the
	matrix). The problem is estimate the number of such singular values to keep. This function uses the theory of Gavish & Donoho (2014) 
	that is valid for matrices with sizes larger than its rank, i.e., for low-rank matrices
	
	Input:
	matrix: matrix to estimate the rank
	
	Output:
	thrKnownNoise: in case sigma is known, multiply this value with sigma and this gives the threshold
	thrUnknownNoise: returns the threshold directly
	noiseEstimation: returns the estimation of the noise
	sv: returns the singular values of the matrix
	"""
	
	m, n = matrix.shape
	beta = 1.0 * m / n
	
	w = (8.0 * beta) / (beta + 1 + np.sqrt(beta**2 + 14 * beta +1))
	lambdaStar = np.sqrt(2.0 * (beta + 1) + w)
	
	omega = 0.56*beta**3 - 0.95*beta**2 + 1.82*beta + 1.43		
	uSVD, wSVD, vSVD = np.linalg.svd(matrix)
	medianSV = np.median(wSVD)
	
	thrKnownNoise = lambdaStar * np.sqrt(n)
	thrUnknownNoise = omega * medianSV
	
	muSqrt = lambdaStar / omega
	noiseEstimation = medianSV / (np.sqrt(n) * muSqrt)		
			
	return thrKnownNoise, thrUnknownNoise, noiseEstimation, wSVD

if __name__ == "__main__":
# Generate a matrix of given rank
	rank = 3
	m = 80
	n = 80
	sigma = 0.02
	X1 = np.random.randn(m,rank)
	X2 = np.random.randn(rank,n)

	Y = np.dot(X1, X2) + sigma * np.random.randn(m,n)

	thrKnownNoise, thrUnknownNoise, noiseEstimation, wSVD = optimalSVHT(Y)
	rankKnownNoise = len(np.where(wSVD > thrKnownNoise*sigma)[0])
	rankUnknownNoise = len(np.where(wSVD > thrUnknownNoise)[0])

	print("Known sigma : {0} - rank={1}".format(thrKnownNoise * sigma, rankKnownNoise))
	print("Unknown sigma : {0} - rank={1} - sigmaNoise={2}".format(thrUnknownNoise, rankUnknownNoise, noiseEstimation))
