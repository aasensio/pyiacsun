from __future__ import print_function
import numpy as np
import scipy.sparse.linalg as splinalg
import scipy.sparse as sp

def alpsLowRank(y, A, At, size, k, tol=1e-3, maxIter=500, frequencyVerbose=20):
	"""
	Recover the low-rank matrix L such that y=A(L)+error
	
	Args:
	    y (real): measurements
	    A (function): forward operator A
	    At (function): adjoint operator At
	    size (int tuple): tuple defining the size of the matrix to recover (m,n)
	    k (int): rank of the matrix to search for
	
	Returns:
	    real array: estimation of the low-rank matrix	
	"""
	yLength = len(y)
	m, n = size
	
	L_cur = np.zeros((m,n))
	QL_cur = np.zeros((m,n))
	QM_cur = np.zeros((m,n))
	
	ALM_cur = np.zeros((yLength,1))
	
	M_cur = np.zeros((m,n))
	M_prev = np.zeros((m,n))
	complementary_Mi = np.ones((m,n))

	I = np.eye(m)
	
	i = 1
	
	while (i <= maxIter):		
		res =  y - A(QL_cur + QM_cur)
		grad = At(res, m, n)
		
## Low rank matrix part
# Active subspace expansion step
		if (i == 1):			
			mat = sp.csr_matrix(grad)
			Uout, _, _ = splinalg.svds(mat, k=k, which='LM', tol=tol)
		else:
			mat = sp.csr_matrix(ortho_UQL_i.dot(grad))
			Uout, _, _ = splinalg.svds(mat, k=k, which='LM', tol=tol)
					
		if (i == 1):
			Si_L = Uout
		else:
			Si_L = np.hstack([Ucur, Uout])			

# Error norm reduction via gradient descent
		proj_grad = Si_L.dot(Si_L.T.dot(grad))		
		mu = np.linalg.norm(proj_grad, 'fro')**2 / np.linalg.norm(A(proj_grad), 2)**2	
		Vi_L = QL_cur + mu*proj_grad

# Best rank-k subspace selection
		mat = sp.csr_matrix(Vi_L)
		UWi_L, SWi_L, VWi_Lt = splinalg.svds(mat, k=k, which='LM')
		SWi_L = np.diag(SWi_L)
		Wi_L = UWi_L.dot(SWi_L.dot(VWi_Lt))
		
# Debias via gradient descene
		res = y-A(Wi_L) - A(M_cur)
		grad = At(res, m, n)
		
		proj_grad = UWi_L.dot(UWi_L.T.dot(grad))
		xi = np.linalg.norm(proj_grad, 'fro')**2 / np.linalg.norm(A(proj_grad), 2)**2
		L_prev = L_cur.copy()
		L_cur = Wi_L + xi*proj_grad
		if (i == 1):
			ULcur = UWi_L.copy()
		else:			
			ULprev = ULcur.copy()
			ULcur = UWi_L.copy()
			
		if (i == 1):
			Ucur = ULcur
		else:
			Ucur = np.hstack([ULcur, ULprev])
			
		ortho_UQL_i = I - Ucur.dot(Ucur.T)

# Update current estimates				
		ALM_prev = np.copy(ALM_cur)
		ALM_cur = A(L_cur + M_prev)
		ALM = ALM_cur - ALM_prev
		tau_LM = ((y-ALM_cur).T.dot(ALM)) / np.linalg.norm(ALM, 2)**2

		QL_cur = L_cur + tau_LM * (L_cur - L_prev)
				
## Sparse matrix part
		# stop()

# Stopping criterion

		norm_cur = np.linalg.norm( (L_cur + M_cur), 'fro')
		norm = np.linalg.norm( (L_cur + M_cur) - (L_prev + M_prev), 'fro')
		normResidual = np.linalg.norm(y - A(L_cur + M_cur))
		if (np.mod(i,frequencyVerbose) == 0):
			print("Iteration {0:4d} : |L(i) - L(i-1)|_F={1:10.7e} - |L(i) - L(i-1)|_F / |L(i-i)|_F={2:10.7e} - residual={3:10.7e}".format(i, \
				norm, norm / norm_cur, normResidual))
		if ( (i > 1) & (norm < tol*norm_cur) ):
			return L_cur

		i += 1
	return L_cur