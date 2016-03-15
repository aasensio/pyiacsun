__all__ = ["prox_rank1_l0"]

from .prox_rank1_generic import prox_rank1_generic
import numpy as np

def hardThreshold(x,t):
	xOut = x.copy()
	xOut[x <= t] = 0.0
	return xOut
	
def prox_rank1_l0(*args, **kwargs):
	"""
	PROX_RANK1_L0 returns the scaled proximity operator for the l0 norm

	x = prox_rank1_l0( x0, D, u )
		where 
		x = argmin_{x} h(x) + 1/2||x-x0||^2_{V}
        and
        V^{-1} = D + u*u'  (or diag(D) + u*u' if D is a vector)
        "D" must be diagonal and positive. "u" can be any vector.
		Here, h(x) = ||x||_0 (the "l-0" norm)
	
	There are also variants:
	x = prox_rank1_l0( x0, D, u, lambda, linTerm, sigma, inverse)
			returns
			x = argmin_{x} h(lambda.*x) + 1/2||x-x0||^2_{V} + linTerm'*x
	        and
			either V^{-1} = D + sigma*u*u' if "inverse" is true (default)
			or     V      = D + sigma*u*u' if "inverse" is false
			and in both cases, "sigma" is either +1 (default) or -1.
			"lambda" should be non-zero
	
	Python : A. Asensio Ramos (March 18, 2015)
	"""
	
	prox = lambda x,t: hardThreshold(x,t)
	
	prox_brk_pts = lambda s : np.hstack((-s,s))
		
	return prox_rank1_generic(prox, prox_brk_pts, *args, **kwargs)