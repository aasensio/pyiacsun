from __future__ import print_function
import numpy as np
from .lte import *
import copy
from ipdb import set_trace as stop

def initLTE(atmos, lines, wavelengthAxis):
    """
    Initialize the LTE synthesis module using nodes
    
    Args:
        atmos (float): array of size (ndepth x 7) defining the reference atmosphere. The columns are
            log(tau)   T [K]   vmic [km/s]   vmac [km/s]   B [G]    thetaB [deg]    phiB [deg]

        lines (float): array of size (nlines x 11) defining the information for the spectral lines. The columns are
            lambda0 [A]    Element    ionization state    log(gf)     Elow [cm^-1]    Lande_up    Lande_low   Jup   Jlow   sigmaABO   alphaABO
            
        wavelengthAxis (float): array of length nlambda that sets the wavelength axis

    Returns:
        None
    """
    initAtmos(atmos)
    initLines(lines, wavelengthAxis)

def synthLTE(referenceAtmos, variablesRF=None, responseFunction=False, deltaRT=0.01):
    """
    Synthesize the Stokes profiles perturbing the reference atmosphere using nodes
    
    Args:
        referenceAtmos (float): array of size (ndepth x 7) defining the reference atmosphere. The columns are
            log(tau)   T [K]   vmic [km/s]   vmac [km/s]   B [G]    thetaB [deg]    phiB [deg]
        variablesRF (optional, list): a list containing (0/1) indicating the variables for which the response functions are obtained
        responseFunction (bool, optional): return the response functions
        deltaRT (float, optional): variation of the parameters when computing the response functions
    
    Returns:
        float: Stokes parameters [4 x nwavelength]
        float: continuum value [nwavelength]        
    """

    logTau = referenceAtmos[:,0]
    
    stokes, cont = synthLines(referenceAtmos)

# Compute the response functions if needed
    if (responseFunction):

        nDepth = len(logTau)
        nLambda = len(cont)    

        atmosPerturbed = np.copy(referenceAtmos)

        typicalValues = [500.0, 1.0, 1.0, 200.0, 50.0, 50.0]

        if (variablesRF == None):
            variablesRF = [1] * 6
        
        nVariables = np.sum(variablesRF)

        RF = np.zeros((nVariables,nDepth,4,nLambda))
        loop = 0
        for indexPar in range(6):            
            if (variablesRF[indexPar] == 1):
                atmosPerturbed = np.copy(referenceAtmos)
                for i in range(nDepth):
                    delta = deltaRT * referenceAtmos[i,indexPar+1]
                    atmosPerturbed[i,indexPar+1] = referenceAtmos[i,indexPar+1] + delta

                    stokesNew, cont = synthLines(atmosPerturbed)

                    atmosPerturbed[i,indexPar+1] = referenceAtmos[i,indexPar+1]

                    RF[loop,i,:,:] = (stokesNew - stokes) / delta

                loop += 1

        return stokes, cont, RF

    else:
        return stokes, cont