from numpy cimport ndarray as ar
from numpy import empty, ascontiguousarray

cdef extern:
	void c_initatmosphere(int *nDepths)
	void c_initlines(int *nLines, double *lineListIn, int *nLambda, double *lambdaIn)
	void c_synthlines(int *nDepths, double *atmosphereIn, int *nLambda, double *stokesOut)
	
cdef int nLambdaGlobal = 1

def initAtmos(ar[double,ndim=2] atmosphereIn):
	cdef:
		int nDepths = atmosphereIn.shape[0]
		ar[double, ndim=2, mode="c"] atmosphere
			
	# Make sure that the 2D array is C_CONTIGUOUS
	atmosphere = ascontiguousarray(atmosphereIn)
	#c_initatmosphere(&nDepths, &atmosphere[0,0])
	c_initatmosphere(&nDepths)
		
def initLines(ar[double,ndim=2] linesIn, ar[double] lambdaIn):
	cdef:
		int nLines = linesIn.shape[0]
		int nLambda = lambdaIn.shape[0]
		ar[double, ndim=2, mode="c"] lines
	global nLambdaGlobal
		
	# Make sure that the 2D array is C_CONTIGUOUS
	lines = ascontiguousarray(linesIn)
	nLambdaGlobal = nLambda
	c_initlines(&nLines, &lines[0,0], &nLambda, &lambdaIn[0])
	
	
def synthLines(ar[double,ndim=2] atmosphereIn):
	global nLambdaGlobal
	cdef:
		int nDepths = atmosphereIn.shape[0]
		int nLambda = nLambdaGlobal
		ar[double, ndim=2, mode="c"] atmosphere
		ar[double, ndim=2] stokesOut = empty((5,nLambdaGlobal), order='F')	
	
	# Make sure that the 2D array is C_CONTIGUOUS
	atmosphere = ascontiguousarray(atmosphereIn)	
	c_synthlines(&nDepths, &atmosphere[0,0], &nLambda, <double*> stokesOut.data)
	
	return stokesOut
