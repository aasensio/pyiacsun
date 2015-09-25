from numpy cimport ndarray as ar
from numpy import empty, ascontiguousarray
import numpy as np

cdef extern:
	void c_setline(int *nLambda, double *lineData, double *waveOut)
	void c_addline(double *lineData)
	void c_milnesynth(int *nLambda, double *modelIn, double *muIn, double *stokesOut)
	void c_milnesynthmany(int *nLambda, int *nModels, double *modelIn, double *muIn, double *stokesOut)

def init(int nLambda, ar[double,ndim=1] lineData):

	cdef:
		ar[double,ndim=1] waveOut = empty(nLambda, order='F')
		
	c_setline(&nLambda, &lineData[0], <double*> waveOut.data)
	
	return waveOut

def addLine(ar[double,ndim=1] lineData):
		
	c_addline(&lineData[0])
	
	return
	
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
	
class milne:
	"""
	Class that synthesizes spectral lines using Milne Eddington.
	To use it:
	from milne import milne as milne
	line = milne(lineInfo)
	lineInfo = [lambda0, JUp, JLow, gUp, gLow, lambdaStart, lambdaStep, nLambda]
	wavelength, stokes = line.synth(model)
	model = [BField, theta, chi, vmac, damping, B0, B1, doppler, kl]
	"""
	
	def __init__(self, nLambda, lineInfo):
		self.lineInfo = lineInfo
		self.nLambda = nLambda		
		self.wavelength = init(nLambda, lineInfo)
		
	
	def synth(self, model, mu=1.0):
		"""
		Synthesize a spectral line using the Milne-Eddington model
		"""		
		stokes = synth(self.nLambda, model, mu)
		
		return stokes
	
	def synthGroup(self, model, mu=1.0):
		"""
		Synthesize many spectral line using the Milne-Eddington model
		"""
		nModels = model.shape[1]		
		stokes = synthGroup(self.nLambda, nModels, model, mu)
		
		return stokes
	
	def __perturbParameter(self, model, index, relativeChange):
		
		newModel = np.copy(model)
		if (model[index] == 0):
			change = relativeChange
		else:
			change = model[index] * relativeChange
		
		newModel[index] += change		
				
		return newModel, change
	
	def __perturbManyParameters(self, model, index, relativeChange):
				
		newModel = np.array(model)		
		change = relativeChange * np.ones(model.shape[1])
		
		ind = np.nonzero(model[index,:])
		change[ind] = newModel[index,ind] * relativeChange
		
		newModel[index,:] += change
				
		return newModel, change
	
	def varChange(self, p):
		model = p
		model[6] = p[6] - p[5]
		model[8] = (p[8] / p[7])**2		
		
		return model
	
	def varInvChange(self, p):
		model = p
		model[6] = p[5] + p[6]		
		model[8] = p[7] * np.sqrt(p[8])
		
		return model
			
	def synthDerivatives(self, model, mu=1.0, relativeChange=1e-3):
		"""
		Compute the derivative of the Stokes profiles with respect to all the variables
		"""
		
		stokes = synth(self.nLambda, model, mu)
		
		stokesDeriv = np.zeros((9,4,self.nLambda))
				
		for i in range(9):			
			newModel, change = self.__perturbParameter(model, i, relativeChange)			
			stokesNew = synth(self.nLambda, newModel, mu)
		
			stokesDeriv[i,:,:] = (stokesNew - stokes) / change			
		
		return stokes, stokesDeriv

	def synthGroupDerivatives(self, model, mu=1.0, relativeChange=1e-3):
		"""
		Compute the derivative of the Stokes profiles with respect to all the variables
		"""
		nModels = model.shape[1]		
		
		stokes = synthGroup(self.nLambda, nModels, model, mu)
				
		stokesDeriv = np.zeros((9,4,self.nLambda,nModels))				
		
		for i in range(9):
			newModel, change = self.__perturbManyParameters(model, i, relativeChange)
						
			stokesNew = synthGroup(self.nLambda, nModels, newModel, mu)
		
# Automatic broadcast "change"
			stokesDeriv[i,:,:,:] = (stokesNew - stokes) / change
		
		return stokes, stokesDeriv
