from numpy cimport ndarray as ar
from numpy import empty, ascontiguousarray

cdef extern from "pymilne.h":
	void c_setline(int *nLambda, double *lineData, double *waveOut)
	void c_milnesynth(int *nLambda, double *modelIn, double *muIn, double *stokesOut)
	void c_milnesynthmany(int *nLambda, int *nModels, double *modelIn, double *muIn, double *stokesOut)

def init(int nLambda, ar[double,ndim=1] lineData):

	cdef:
		ar[double,ndim=1] waveOut = empty(nLambda, order='F')
		
	c_setline(&nLambda, &lineData[0], <double*> waveOut.data)
	
	return waveOut
	
def synth(int nLambda, ar[double,ndim=1] modelIn, double muIn):
	
	cdef:
		ar[double,ndim=2] stokesOut = empty((4,nLambda), order='F')
		
	c_milnesynth(&nLambda, &modelIn[0], &muIn, <double*> stokesOut.data)
	
	return stokesOut
	
def synthGroup(int nLambda, int nModels, ar[double,ndim=2] modelIn, double muIn):
	
	cdef:
		ar[double,ndim=3] stokesOut = empty((4,nLambda,nModels), order='F')
		ar[double, ndim=2, mode="c"] model
		
	# Make sure that the 2D array is C_CONTIGUOUS
	model = ascontiguousarray(modelIn)
		
	c_milnesynthmany(&nLambda, &nModels, &model[0,0], &muIn, <double*> stokesOut.data)
	
	return stokesOut