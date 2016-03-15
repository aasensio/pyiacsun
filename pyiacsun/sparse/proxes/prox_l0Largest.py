import numpy as np

def prox_l0Largest(x, lambdaPar):
    """
    Hard thresholding operator
    """
    xPar = np.copy(x)
    ind = np.argsort(np.abs(x))
    xPar[ind[lambdaPar:]] = 0.0
    return xPar