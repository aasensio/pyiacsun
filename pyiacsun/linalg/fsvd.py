__all__ = ['fsvd']
import numpy as np

def fsvd(A, k, i, usePowerMethod=False):

#function [U,S,V] = fsvd(A, k, i, usePowerMethod)
# FSVD Fast Singular Value Decomposition 
# 
#   [U,S,V] = FSVD(A,k,i,usePowerMethod) computes the truncated singular
#   value decomposition of the input matrix A upto rank k using i levels of
#   Krylov method as given in [1], p. 3.
# 
#   If usePowerMethod is given as true, then only exponent i is used (i.e.
#   as power method). See [2] p.9, Randomized PCA algorithm for details.
# 
#   [1] Halko, N., Martinsson, P. G., Shkolnisky, Y., & Tygert, M. (2010).
#   An algorithm for the principal component analysis of large data sets.
#   Arxiv preprint arXiv:1007.5510, 0526. Retrieved April 1, 2011, from
#   http://arxiv.org/abs/1007.5510. 
#   
#   [2] Halko, N., Martinsson, P. G., & Tropp, J. A. (2009). Finding
#   structure with randomness: Probabilistic algorithms for constructing
#   approximate matrix decompositions. Arxiv preprint arXiv:0909.4061.
#   Retrieved April 1, 2011, from http://arxiv.org/abs/0909.4061.
# 
#   See also SVD.
# 
#   Copyright 2011 Ismail Ari, http://ismailari.com.

	if (usePowerMethod == False):
		i = 1;
		
	s = A.shape

# Take (conjugate) transpose if necessary. It makes H smaller thus
# leading the computations to be faster
	isTransposed = False
	if (s[0] < s[1]):
		A = A.T
		isTransposed = True

	n = A.shape[1]
	l = k + 2

# Form a real nxl matrix G whose entries are iid Gaussian r.v.s of zero
# mean and unit variance
	G = np.random.randn(n,l)		

	if (usePowerMethod):

# Use only the given exponent
		H = np.dot(A,G)
		for j in range(2,i+1):
			H = np.dot(A, np.dot(A.T,H))
	else:
		
# Compute the mxl matrices H^{(0)}, ..., H^{(i)}
# Note that this is done implicitly in each iteration below.
		H = []
		H = np.append(A*G)
		for j in range(1,i):
			H = np.append(np.dot(A, np.dot(A.T, H[j-1])))
		H = np.concatenate(H)
			
## Using the pivoted QR-decomposiion, form a real mx((i+1)l) matrix Q
## whose columns are orthonormal, s.t. there exists a real
## ((i+1)l)x((i+1)l) matrix R for which H = QR.      
	[Q, R] = np.linalg.qr(H)
	
	#pdb.set_trace()

## Compute the nx((i+1)l) product matrix T = A^T Q
	T = np.dot(A.T,Q)

## Form an SVD of T
	Vt, St, W = np.linalg.svd(T)

## Compute the mx((i+1)l) product matrix
	Ut = np.dot(Q,W)

## Retrieve the leftmost mxk block U of Ut, the leftmost nxk block V of
## Vt, and the leftmost uppermost kxk block S of St. The product U S V^T
## then approxiamtes A. 

	if (isTransposed):
		V = Ut[:,0:k-1];
		U = Vt[:,0:k-1];     
	else:
		U = Ut[:,0:k-1]
		V = Vt[:,0:k-1]
	
	S = St[0:k-1]
	
	return U, S, V

#A = np.array([[1,2,3,4],[3,5,6,7],[2,8,3,1]])
#A = np.random.randn(2000,2000)

#start = time.clock()
#U2, S2, V2 = np.linalg.svd(A)
#print (time.clock() - start)

#start = time.clock()
#U, S, V = fsvd(A, 9, 3, usePowerMethod=True)
#print (time.clock() - start)
