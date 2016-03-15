__all__ = ["prox_rank1_box"]

from .prox_rank1_generic import prox_rank1_generic
import numpy as np

def prox_rank1_box(lwr, upr, *args, **kwargs):
	"""
	PROJ_RANK1_BOX returns the scaled proximity operator for box constraints

	x = proj_rank1_box( lwr, upr, x0, D, u )
			where 
			x = argmin_{x} h(x) + 1/2||x-x0||^2_{V}
			and
			V^{-1} = D + u*u'  (or diag(D) + u*u' if D is a vector)
			"D" must be diagonal and positive. "u" can be any vector.

	Here, h(x) is the indicator function of the set
		{ x : lwr <= x <= upr }
		(Set any component of lwr to -Inf and upr to +Inf to effectively
			ignore those particular constraints)

	There are also variants:
	x = proj_rank1_box( lwr, upr, x0, D, u, lambda, linTerm, sigma, inverse)
			returns
			x = argmin_{x} h(lambda.*x) + 1/2||x-x0||^2_{V} + linTerm'*x
			and
			either V^{-1} = D + sigma*u*u' if "inverse" is true (default)
			or     V      = D + sigma*u*u' if "inverse" is false
			and in both cases, "sigma" is either +1 (default) or -1.
			"lambda" should be non-zero

	Note that UNLIKE prox_rank1_l1.m and other functions, the calling
	sequence is slighty different, since you must pass in "lwr" and "upr"

	Stephen Becker, Feb 26 2014, stephen.beckr@gmail.com
	Reference: "A quasi-Newton proximal splitting method" by S. Becker and J. Fadili
	NIPS 2012, http://arxiv.org/abs/1206.1156

	BSee also prox_rank1_generic.m
	
	Python version: A. Asensio Ramos (March 12, 2015)
	"""
	
	prox = lambda x, t : np.max([np.min([upr,x]), lwr])
	prox_brk_pts = lambda s : np.hstack((lwr,upr))
		
	return prox_rank1_generic(prox, prox_brk_pts, *args, **kwargs)
