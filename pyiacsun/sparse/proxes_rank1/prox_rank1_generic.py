__all__ = ["prox_rank1_generic"]

import numpy as np
import sys

def prox_rank1_generic(prox, prox_brk_pts, x0, d, u = None, Lambda = None, linTerm = None, plusminus=1, invert = True, verbose = False):
	"""
	PROX_RANK1_GENERIC returns the scaled proximity operator for a generic function h
	(provided the generic function is separable and has a piece-wise linear prox)
	This function is intended be used as follows:

	(1) Instantiate:
		scaledProx = lambda x0, d, u, varargin = None : prox_rank1_l1(x0, d, u, l)@(varargin) prox_rank1_generic( prox, prox_brk_pts,varargin{:})
			where 'prox' and 'prox_brk_pts' implicitly define the function h
			i.e., prox(x0,t) = argmin_{x} t*h(x) + 1/2||x-x0||^2
			and 
				prox_brk_pts(t) returns a row-vector with the break points
				that specify where t*h(x) is piecewise linear
	            (this is if h(x) = [ h_1(x_1); ... ; h_n(x_n) ]. If instead not
					all the h_i are identical, prox_brk_pts(t) should return
					a matrix).
			See the examples below because prox_brk_pts must allow a vector "t"
				so you must define this appropriately.
	
	(2) Call the "scaledProx" function, which has signature:
		x = scaledProx( x0, D, u )
			where 
			x = argmin_{x} h(x) + 1/2||x-x0||^2_{V}
			and
			V^{-1} = D + u*u'  (or diag(D) + u*u' if D is a vector)
			"D" must be diagonal and positive. "u" can be any vector.

		There are also variants:
	
		x = scaledProx( x0, D, u, lambda, linTerm, sigma, inverse)
			returns
			x = argmin_{x} h(lambda.*x) + 1/2||x-x0||^2_{V} + linTerm'*x
			and
			either V^{-1} = D + sigma*u*u' if "inverse" is true (default)
			or     V      = D + sigma*u*u' if "inverse" is false
			and in both cases, "sigma" is either +1 (default) or -1.
			"lambda" should be non-zero
	
	Examples:
      1. if h(x) = ||x||_1 then
           prox            = @(x,t) sign(x).*max(0, abs(x) - t );
           prox_brk_pts    = @(t) [-t,t];
      2. if h(x) is the indicator function of the set { x : x >= 0}, then
           prox            = @(x,t) max(0, x);
           prox_brk_pts    = @(t) 0; 
      3. if h(x) is the indicator function of the set { x : lwr <= x <= upr }
           where lwr and upr are vectors, then
           prox            = @(x,t) max( min(upr,x), lwr );
           prox_brk_pts    = @(t) [lwr,upr];  (Note: this is a matrix)
      4. if h(x) is the hinge-loss h(x) = max( 1-x, 0 ), then
           prox        = @(x,t) 1 + (x-1).*( x > 1 ) + (x + t - 1).*( x + t < 1  );
           prox_brk_pts    = @(t)[ones(size(t)), 1-t];
      5. if h(x) is the indicator function of the l_infinity ball, then
           prox            = @(x,t) sign(x).*min( 1, abs(x) );
           prox_brk_pts    = @(t) [-ones(size(t)),ones(size(t))]; 


	Stephen Becker, Feb 26 2014, stephen.beckr@gmail.com
	Reference: "A quasi-Newton proximal splitting method" by S. Becker and J. Fadili
	NIPS 2012, http://arxiv.org/abs/1206.1156
	
	Python version directly translated from Matlab version (including comments): A. Asensio Ramos (March 12, 2015)
	"""

	if (len(u) == 0):
		u = np.asarray([[0]])
		
	if (np.all(u == 0)):
# Diagonal scaling
		NO_U = True
	else:
		NO_U = False

	if (NO_U):
		uinv = 0
	else:
		uinv = u / d / np.sqrt(1.0 + u.T.dot(u/d))

# Check for positive definiteness of matrix V in a special case (warning: not in all cases)
	if ((plusminus < 0) & (np.all(d == d[0]))):
		minE = d[0] + plusminus * np.linalg.norm(u)
		if (minE <= 0):
			print("The scaling matrix is not positive definite")
			sys.exit(1)
		
# In all cases, we find prox_h^V, but how we define V
# in terms of d and u depends on "INVERT"
	if (invert):
# So V^{-1} = diag(d)     + sigma*u*u'
# and     V = diag(1./d)  - sigma*uinv*uinv'
		Vinv = lambda y : d / y + plusminus * (u.T.dot(y)).dot(u)
		dd = 1.0 / d
		uu = uinv
		plusminus = -plusminus
		
		if (NO_U):
			ud = u / np.sqrt(1.0+u*(u/d))
		else:
			ud = u / np.sqrt(1.0+u.T.dot(u/d))
		dInv = 1.0 / dd
		
	else:
# So V = diag(d)     + sigma*u*u'
# and     V^(-1) = diag(1./d)  - sigma*uinv*uinv'
		Vinv = lambda y : y / d - plusminus * (uinv.T.dot(y)).dot(uinv)
		dd = d.copy()
		uu = u.copy()
		
		ud = uu / dd
		dInv = 1.0 / dd

	if (NO_U):
		uu = 0
		ud = 0
	
# We make a change of variables, e.g., x <-- lambda*.x
# change x0 <-- lambda.*x0, linTerm <-- linTerm./lambda
# and V <-- diag(1./lambda)*V*diag(1./lambda). Because V is defined
# implicitly, and depends on INVERT, this is a bit of a headache.
# We'll do some changes here, and some later in the code
# e.g., combine linTerm and V scaling so we don't have to redefine Vinv
	if (Lambda != None):
		if (np.any(Lambda) == 0):
			print("scaling factor lambda must be non-zero")

# note that lambda < 0 should be OK
		x0 *= Lambda

# Scale V = diag(dd) + sigma*uu*uu' by V <-- diag(1./lambda)*V*diag(1./lambda)		
		dd /= Lambda**2
		uu /= Lambda
		ud *= Lambda
		dInv = 1.0 / dd
	
	t = prox_brk_pts(1.0 / dd)
	
	if ((linTerm != None)):
		if (np.linalg.norm(linTerm) >= 0):
			if (Lambda != None):
				x0 -= Vinv(linTerm)
			else:
# V is scaled V <-- diag(1./lambda)*V*diag(1./lambda)
#   so Vinv is scaled the opposite.
# linTerm is scaled linTerm <== linTerm./lambda % V is scaled V <-- diag(1./lambda)*V*diag(1./lambda)
#   so Vinv is scaled the opposite.
# linTerm is scaled linTerm <== linTerm./lambda
				x0 -= Lambda * Vinv(linTerm)		
	
# The main heart
	X = lambda a : prox(x0 - plusminus * a * ud, dInv)
	
# Only return if we have only diagonal scaling
	if (NO_U):		
# In this case, alpha is irrelevant
		x = prox(x0, dInv)
		if (Lambda != None):
# Undo the lambda scaling
			x /= Lambda
		return x
	
	brk_pts = plusminus * (dd/uu) * (x0-t)
	brk_pts = np.unique(brk_pts)
	brk_pts = brk_pts[np.where(np.isfinite(brk_pts))]
	
# Main loop
# Lower bound are a for which p<=0
# Upper bound are a for which p>0
# If a is increasing, so is p(a)
	lwrBnd = 0
	uprBnd = len(brk_pts)
	iMax = int(np.ceil(np.log2(len(brk_pts))) + 1)
	for i in range(iMax):
		if (uprBnd - lwrBnd <= 1):
			if (verbose):
				print("Bounds are too close")
			break
		j = int(round(np.mean([lwrBnd,uprBnd])))
		if (verbose):
			print("j is {0} (bounds were [{1},{2}])".format(j,lwrBnd,uprBnd))
		
		if (j == lwrBnd):
			j += 1
		elif (j == uprBnd):
			j -= 1
		
		a = brk_pts[j]
		x = X(a)
		
		p = a + uu.T.dot(x0-x)
		
		if (p > 0):
			uprBnd = j
		elif (p < 0):
			lwrBnd = j
			
		
	cnt = i  # Number of iterations
	
# Now, determine linear part, which we infer from two points.
# If lwr/upr bounds are infinite, we take special care
# e.g., we make a new "a" slightly lower/bigger, and use this
# to extract linear part.
	if (lwrBnd == 0):
		a2 = brk_pts[uprBnd]
		a1 = a2 - 10     # Arbitrary
		aBounds = np.asarray([-np.inf, a2])
	elif (uprBnd == len(brk_pts)):
		a1 = brk_pts[lwrBnd]
		a2 = a1 + 10     # Arbitrary
		aBounds = np.asarray([a1,np.inf])
	else:
# In general case, we can infer linear part from the two break points
		a1 = brk_pts[lwrBnd]
		a2 = brk_pts[uprBnd]
		aBounds = np.asarray([a1, a2])
		
	x1 = X(a1)
	x2 = X(a2)
	dx = (x2-x1) / (a2-a1)
	
# Thus for a in (a1,a2), x(a) = x1 + (a-a1)*dx
# Solve 0 = a + dot( uu, y - (x1 + (a-a1)*dx ) )
#         = a + dot(uu,y - x1 + a1*dx ) - a*dot(uu,dx)
# so:
	a = uu.T.dot(x0 - x1 + a1 * dx) / (-1.0 + uu.T.dot(dx))
	
	if ((a < aBounds[0]) | (a > aBounds[1])):
		print("alpha is not the correct range")
		
# The solution
	x = X(a)
	
	if (Lambda != None):
		x /= Lambda
	
	if (verbose):
		print("Took {0} of {1} iterations, lwrBnd is {2} {3}".format(i,iMax,lwrBnd, len(brk_pts)))
		
	return x