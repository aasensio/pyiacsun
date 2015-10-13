__all__ = ['rebin']
import numpy as np
# http://wiki.scipy.org/Cookbook/Rebinning
def rebin(a, *args):
	'''rebin ndarray data into a smaller ndarray of the same rank whose dimensions
	are factors of the original dimensions. eg. An array with 6 columns and 4 rows
	can be reduced to have 6,3,2 or 1 columns and 4,2 or 1 rows.
	
	Args:
	    a (array): array to be rebinned
	    *args: dimensions to rebin the array to
	
	example usages:
	>>> a=rand(6,4); b=rebin(a,3,2)
	>>> a=rand(6); b=rebin(a,2)
	'''
	shape = a.shape
	lenShape = len(shape)
	factor = np.asarray(shape)/np.asarray(args)
	evList = ['a.reshape('] + \
             ['args[%d],factor[%d],'%(i,i) for i in range(lenShape)] + \
             [')'] + ['.sum(%d)'%(i+1) for i in range(lenShape)] + \
             ['/factor[%d]'%i for i in range(lenShape)]	
	return eval(''.join(evList))

