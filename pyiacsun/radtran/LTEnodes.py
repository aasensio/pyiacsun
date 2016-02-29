from __future__ import print_function
import numpy as np
from .lte import *
import copy

__all__ = ['initLTENodes', 'nodePositions', 'synthLTENodes']
def initLTENodes(atmos, lines, wavelengthAxis):
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

def _interpolateNodes(logTau, variable, nodes):
    n = logTau.shape[0]
    out = np.zeros_like(variable)

    if (len(nodes) == 0):
        out = variable

    if (len(nodes) == 1):
        out = variable + nodes[0]
    
    if (len(nodes) >= 2):
        pos = np.linspace(0, n-1, len(nodes), dtype=int)
        coeff = np.polyfit(logTau[pos], nodes, len(nodes)-1)
        out = variable + np.polyval(coeff, logTau)
    return out

def nodePositions(logTau, nodes):
    nDepth = logTau.shape[0]
    positions = []

    for n in nodes:
        if (len(n) == 0):
            positions.append([])

        if (len(n) == 1):
            positions.append([int(nDepth/2)])
    
        if (len(n) >= 2):
            pos = np.linspace(0, nDepth-1, len(n), dtype=int)
            positions.append(pos)

    return positions


def synthLTENodes(referenceAtmos, nodes, responseFunction=False, deltaRT=0.01):
    """
    Synthesize the Stokes profiles perturbing the reference atmosphere using nodes
    
    Args:
        referenceAtmos (float): array of size (ndepth x 7) defining the reference atmosphere. The columns are
            log(tau)   T [K]   vmic [km/s]   vmac [km/s]   B [G]    thetaB [deg]    phiB [deg]
        nodes (list): a list containing six sets of nodes for each physical variable. Depending on the length of each
            sublist, the perturbation will be applied using different functional forms: 
                - constant if only one node is used
                - linear perturbing the values at the extremes if two nodes are used
                - polynomial at equispaced positions in the log(tau) axis if more than 2 nodes are used
        responseFunction (bool, optional): return the response functions
        deltaRT (float, optional): variation of the parameters when computing the response functions
    
    Returns:
        float: Stokes parameters [4 x nwavelength]
        float: continuum value [nwavelength]
        float: modified atmosphere [ndepth x 7]
    """

    logTau = referenceAtmos[:,0]
    atmos = np.copy(referenceAtmos)
    
    for i, n in enumerate(nodes):
        atmos[:,i+1] = _interpolateNodes(logTau, referenceAtmos[:,i+1], n)

    stokes, cont = synthLines(atmos)

    atmosPerturbed = np.copy(referenceAtmos)

    typicalValues = [500.0, 1.0, 1.0, 200.0, 50.0, 50.0]

    RF = []

# Compute the response functions if needed
    if (responseFunction):
        for indexPar in range(6):
            RFNode = []
            for indexNode in range(len(nodes[indexPar])):
                nodesNew = copy.deepcopy(nodes)
                nodesNew[indexPar][indexNode] += deltaRT * typicalValues[indexPar]

                delta = nodesNew[indexPar][indexNode] - nodes[indexPar][indexNode]

                for i, n in enumerate(nodesNew):
                    atmosPerturbed[:,i+1] = _interpolateNodes(logTau, referenceAtmos[:,i+1], n)

                stokesNew, cont = synthLines(atmosPerturbed)

                RFNode.append((stokesNew - stokes) / delta)

            RF.append(RFNode)
        return stokes, cont, atmos, RF

    else:
        return stokes, cont, atmos