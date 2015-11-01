__all__ = ['lineMin']

import numpy as np
import scipy.optimize as op

def lineMin(spectrum, x=None, findMax=False, estimated=None, rangeAround=None, rangeFit=3, order=4):
    """
    Find the minimum of a spectral line with subpixel precision
    
    Args:
        spectrum (TYPE): spectrum containing spectral lines
        x (float, optional): x-axis
        findMax (bool, optional): by default, it finds the minimum. If this is set, it detects the maximum
        estimated (TYPE, optional): first estimate of the position index in the spectrum array where (default: size/2)
        rangeAround (TYPE, optional): range around the estimated position where the minimum is seeked (default: size/4)
        rangeFit (int, optional): once a crude estimation is found, range around this estimation that is used for the refinement (default: 3)
        order (int, optional): order of the polynomial used for refinement (default: 4)
    
    Returns:
        TYPE: Description
    
    Deleted Args:
        x (TYPE, optional): Description
    """

    nx = spectrum.shape[0]
    if (estimated == None):
        estimated = nx / 2
    if (rangeAround == None):
        rangeAround = estimated
    xAxis = np.arange(nx)

    low = np.max([0, int(estimated - rangeAround)])
    up = np.min([int(estimated + rangeAround), nx])

    if (findMax):
        loc = np.argmax(spectrum[low:up]) + low
        factor = -1.0        
    else:
        loc = np.argmin(spectrum[low:up]) + low
        factor = 1.0

    coeff = np.polyfit(xAxis[loc-rangeFit:loc+rangeFit], spectrum[loc-rangeFit:loc+rangeFit] * factor, order)
    p = np.poly1d(coeff)

    out = op.minimize_scalar(p, method='Bounded', bounds=(xAxis[loc-rangeFit], xAxis[loc+rangeFit])).x    
    
    if (x != None):
        out = np.interp(out, xAxis, x)

    return out


# if (__name__ == '__main__'):
#     x = np.linspace(0.0,30.0,150)
#     y = 1.0 - np.exp(-(x-15.0)**2)

#     res = lineMin(y, x=x, estimated=60)

#     # pl.plot(x,y)
#     # pl.axvline(res.x)
#     print(res)
