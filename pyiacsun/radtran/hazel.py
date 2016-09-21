import numpy as np
from .pyhazel import *

__all__ = ["hazel"]

class hazel:
    """
    Class that synthesizes spectral lines using Hazel
    To use it:

    out = hazel()
    l, stokes, dStokes = hazel.synthDerivatives(whichPars, 1e-3, *args)
    with *args a tuple containing all parameters described in the synth method

    """

    def __init__(self):
        c_init()

    def synth(self, *args):     
        """
        Carry out a synthesis with Hazel (see the manual for the meaning of all of them)
        
        Args: 
            synModeInput: (int) synthesis mode
            nSlabsInput: (int) number of slabs
            B1Input: (float) vector of size 3 with the magnetic field vector in spherical coordinates for the first component
            B2Input: (float) vector of size 3 with the magnetic field vector in spherical coordinates for the second component
            hInput: (float) height
            tau1Input: (float) optical depth of the first component
            tau2Input: (float) optical depth of the second component        
            boundaryInput: (float) vector of size 4 with the boundary condition for (I,Q,U,V)
            transInput: (int) transition to compute from the model atom
            atomicPolInput: (int) include or not atomic polarization
            anglesInput: (float) vector of size 3 describing the LOS
            nLambdaInput: (int) number of wavelength points
            lambdaAxisInput: (float) vector of size nLambdaInput with the wavelength axis (relative to 10829.0911 A)        
            dopplerWidth1Input: (float) Doppler width of the first component
            dopplerWidth2Input: (float) Doppler width of the second component
            dampingInput: (float) damping
            dopplerVelocityInput: (float) bulk velocity affecting the first component
            dopplerVelocity2Input: (float) bulk velocity affecting the second component
            ffInput: (float) filling factor
            betaInput: (float) value to be multiplied by the source function of the second component to allow for emission lines in the disk
            nbarInput: (float) vector of size 4 to define nbar for every transition of the model atom (set them to zero to use Allen's)
            omegaInput: (float) vector of size 4 to define omega for every transition of the model atom (set them to zero to use Allen's)
            
        Returns:
            wavelengthOutput: (float) vector of size nLambdaInput with the wavelength axis
            stokesOutput: (float) array of size (4,nLambdaInput) with the emergent Stokes profiles
            epsOutput: (float) array of size (4,nLambdaInput) with the emissivity vector at each wavelength
            etaOutput: (float) array of size (4,4,nLambdaInput) with the propagation matrix at each wavelength
        """             
        return _synth(*args)

    def __synthPerturbation(self, paramToModify, perturbation, *args):

        newPars = list(args)

        if (paramToModify == 'B1'):
            newPars[2][0] += perturbation
        if (paramToModify == 'thB1'):
            newPars[2][1] += perturbation
        if (paramToModify == 'chiB1'):
            newPars[2][2] += perturbation

        if (paramToModify == 'B2'):
            newPars[3][0] += perturbation
        if (paramToModify == 'thB2'):
            newPars[3][1] += perturbation
        if (paramToModify == 'chiB2'):
            newPars[3][2] += perturbation

        if (paramToModify == 'tau1'):
            newPars[5] += perturbation
        if (paramToModify == 'tau2'):
            newPars[6] += perturbation

        if (paramToModify == 'vth1'):
            newPars[13] += perturbation
        if (paramToModify == 'vth2'):
            newPars[14] += perturbation

        if (paramToModify == 'a'):
            newPars[15] += perturbation

        if (paramToModify == 'v1'):
            newPars[16] += perturbation
        if (paramToModify == 'v2'):
            newPars[17] += perturbation

        if (paramToModify == 'ff'):
            newPars[18] += perturbation

        if (paramToModify == 'beta'):
            newPars[19] += perturbation

        newPars = tuple(newPars)

        _, stokesPerturbed, _, _ = _synth(*newPars)

        return stokesPerturbed


    def synthDerivatives(self, paramsRF, perturbation, *args):
        """
        Carry out a synthesis with Hazel and also compute response functions to the parameters
        
        Args: 
            (see the manual for the meaning of all of them)
            - paramsRF : (str) list containing any combination of the following parameters to which the RF are computed
                ['B1', 'thB1', 'chiB1', 'B2', 'thB2', 'chiB2', 'tau1', 'tau2', 'vth1','vth2','a','v1','v2', 'ff', 'beta']
            - perturbation: (float) perturbation used for computing the numerical RF
            - synModeInput: (int) synthesis mode
            - nSlabsInput: (int) number of slabs
            - B1Input: (float) vector of size 3 with the magnetic field vector in spherical coordinates for the first component
            - B2Input: (float) vector of size 3 with the magnetic field vector in spherical coordinates for the second component
            - hInput: (float) height
            - tau1Input: (float) optical depth of the first component
            - tau2Input: (float) optical depth of the second component        
            - boundaryInput: (float) vector of size 4 with the boundary condition for (I,Q,U,V)
            - transInput: (int) transition to compute from the model atom
            - atomicPolInput: (int) include or not atomic polarization
            - anglesInput: (float) vector of size 3 describing the LOS
            - nLambdaInput: (int) number of wavelength points
            - lambdaAxisInput: (float) vector of size nLambdaInput with the wavelength axis (relative to 10829.0911 A)        
            - dopplerWidth1Input: (float) Doppler width of the first component
            - dopplerWidth2Input: (float) Doppler width of the second component
            - dampingInput: (float) damping
            - dopplerVelocityInput: (float) bulk velocity affecting the first component
            - dopplerVelocity2Input: (float) bulk velocity affecting the second component
            - ffInput: (float) filling factor
            - betaInput: (float) value to be multiplied by the source function of the second component to allow for emission lines in the disk
            - nbarInput: (float) vector of size 4 to define nbar for every transition of the model atom (set them to zero to use Allen's)
            - omegaInput: (float) vector of size 4 to define omega for every transition of the model atom (set them to zero to use Allen's)
            
        Returns:
            wavelengthOutput: (float) vector of size nLambdaInput with the wavelength axis
            stokesOutput: (float) array of size (4,nLambdaInput) with the emergent Stokes profiles
            stokesDeriv: (float) array of size (nPar,4,nLambdaInput) with the response functions or each indicated parameter
        """ 

        wavelengthOutput, stokes, epsOutput, etaOutput = _synth(*args)

        nRF = len(paramsRF)

        stokesDeriv = np.zeros((nRF,4,len(wavelengthOutput)))

        for i in range(nRF):
            stokesNew = self.__synthPerturbation(paramsRF[i], perturbation, *args)

            stokesDeriv[i,:,:] = (stokesNew - stokes) / perturbation            

        return wavelengthOutput, stokes, stokesDeriv
