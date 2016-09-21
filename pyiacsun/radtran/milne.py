from __future__ import print_function
import numpy as np
from .pymilne import *

__all__ = ["milne"]

class milne:
    """
    Class that synthesizes spectral lines using a Milne Eddington model atmosphere.
    To use it::
    
        from milne import milne as milne
        
        lineInfo = [lambda0, JUp, JLow, gUp, gLow, lambdaStart, lambdaStep, nLambda]
        line = milne(lineInfo)
        
        model = [BField, theta, chi, vmac, damping, B0, B1, doppler, kl]
        wavelength, stokes = line.synth(model)

    The class needs to be initalized with a spectral line. In case the class is initialized without
    any line information, it will list possible spectral lines. In such case, lines can be later
    added using::

        lineInfo = [lambda0, JUp, JLow, gUp, gLow, lambdaStart, lambdaStep, nLambda]
        out.addLine(lineInfo)
    """
    
    def __init__(self, nLambda=None, lineInfo=None):
        if ((nLambda == None) and (lineInfo == None)):
            print("""
            The MILNE class needs as input the number of wavelengths and the line information:")

            out = milne(nLambda, lineInfo)

            with lineInfo = [lambda0, JUp, JLow, gUp, gLow, lambdaStart, lambdaStep, nLambda]

            A sample of possible lines follows:

            0  6301.5080  1.833  2.0  1.500  2.0  6300.80  6301.99  0.001
            1  6302.4904  2.500  1.0  0.000  0.0  6302.00  6303.19  0.001
            2  5247.0500  1.500  2.0  1.746  3.0  5246.50  5247.50  0.001
            3  5250.2090  0.000  0.0  3.000  1.0  5249.80  5250.70  0.001 
            4  5247.5650  0.000  0.0  2.500  1.0  5247.00  5247.90  0.001
            5  6173.3340  2.500  1.0  0.000  0.0  6173.00  6173.80  0.001
            6  7090.3816  1.000  0.0  3.000  1.0  7089.80  7091.00  0.001

            Afterwards, lines can be added to the problem for blends using

            out.addLine(lineInfo)

            """)
        else:
            self.lineInfo = lineInfo
            self.nLambda = nLambda      
            self.wavelength = init(nLambda, np.asarray(lineInfo))

    def addLine(self, lineInfo):
        """
        Add a spectral line to the list of lines that will be synthesized.
        Usage::

            lineInfo = [lambda0, JUp, JLow, gUp, gLow, lambdaStart, lambdaStep, nLambda]
            milne.addLine(lineInfo)
        """
        addLine(np.asarray(lineInfo))
        return
    
    def synth(self, model, mu=1.0):
        """Synthesize a spectral line using the Milne-Eddington model and the
        model passed as parameter. The model is given by a list containing the 
        following information in order. The source function is linear with optical
        depth, so that S=B0+B1*tau

        Args:
            BField (float): magnetic field strength [G]
            theta (float): inclination of the magnetic field in the line-of-sight reference system [deg]
            chi (float): azimuth of the magnetic field in the line-of-sight reference system  [deg]
            vmac (float): bulk velocity [km/s]
            damping (float): damping parameter of the line
            B0 (float): parameter defining the source function
            B1 (float): parameter defining the source function
            doppler (float): width of the line [milliAngstrom]
            kl: (float): line/continuum opacity ratio

        Output:
            stokes (float): 4 x nLambda array with the Stokes parameters (I,Q,U,V)
        """     
        stokes = synth(self.nLambda, model, mu)
        
        return stokes
    
    def synthGroup(self, model, mu=1.0):
        """Synthesize a spectral line using the Milne-Eddington model like ``synth`` but
        for a set of models. This is faster than calling many times to ``synth``.

        Args:
            model (float): array of size 9 x nModels with all the model parameters as defined
                in ``synth`` for all cases.

        Returns:
            stokes (float): 4 x nLambda x nModels array with the Stokes parameters (I,Q,U,V)
        
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
            
    def synthDerivatives(self, model, mu=1.0, relativeChange=1e-3):
        """Synthesize Stokes parameters and compute the derivative of the Stokes profiles with 
        respect to all the variables 

        Args:
            model (float): array of size 9 with all the model parameters as defined
                in ``synth`` for all cases.
            relativeChange (optional): relative change for all variables to compute the derivatives.
                They are computed numerically. Default value: 1e-3

        Returns:
            stokes (float): 4 x nLambda array with the Stokes parameters (I,Q,U,V)
            stokesDeriv (float): 9 x 4 x nLambda array with the derivative of the Stokes parameters (I,Q,U,V)
                with respect to all model parameters
        
        """        
        stokes = synth(self.nLambda, model, mu)
        
        stokesDeriv = np.zeros((9,4,self.nLambda))
                
        for i in range(9):          
            newModel, change = self.__perturbParameter(model, i, relativeChange)            
            stokesNew = synth(self.nLambda, newModel, mu)
        
            stokesDeriv[i,:,:] = (stokesNew - stokes) / change          
        
        return stokes, stokesDeriv

    def synthGroupDerivatives(self, model, mu=1.0, relativeChange=1e-3):
        """Synthesize Stokes parameters and compute the derivative of the Stokes profiles with 
        respect to all the variables for many models

        Args:
            model (float): array of size 9 x nModels with all the model parameters as defined
                in ``synth`` for all cases.
            relativeChange (optional): relative change for all variables to compute the derivatives.
                They are computed numerically. Default value: 1e-3

        Returns:
            stokes (float): 4 x nLambda x nModels array with the Stokes parameters (I,Q,U,V)
            stokesDeriv (float): 9 x 4 x nLambda x nModels array with the derivative of the Stokes parameters (I,Q,U,V)
                with respect to all model parameters
        
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
