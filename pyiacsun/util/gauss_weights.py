__all__ = ['gauss_weights']

import numpy as np

def gauss_weights(interval, n):
    """
    Return the location and weights for the integration of a function using the Gauss-Legendre function
    
    Args:
        interval (tuple): range of integration
        n (int): order of the integration
    
    Returns:
        x, w: position and weightrs
        
    """
    x, w = np.polynomial.legendre.leggauss(n)

    t = 0.5*(x+1)*(interval[1]-interval[0]) + interval[0]
    w *= 0.5*(interval[1]-interval[0])

    return t, w