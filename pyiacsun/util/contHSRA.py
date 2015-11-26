__all__ = ['contHSRA']

import numpy as np

def contHSRA(wavelength, units='frequency'):
    """
    Return the continuum spectrum of the HSRA model, obtained as a 3rd degree polynomial fit. This gives
    a maximum error of 0.125%
    
    Args:
        wavelength (float): wavelength in Angstrom
        units (str, optional): 'frequency' or 'wavelength' units
    
    Returns:
        float: value of the continuum in frequency units
    """

    spec = np.zeros_like(wavelength)
    c = np.asarray([-4.906054765549e13,1.684734544039e11,1.507254517567e7,-7561.242976546])
    mask = wavelength < 3644.15
    w = wavelength[mask]
    spec[mask] = np.polyval(c[::-1], w)
    
    c = np.asarray([-4.4650822755e14,6.1319780351059e11,-9.350928003805e7])
    mask = (wavelength >= 3644.15) & (wavelength < 3750)
    w = wavelength[mask]
    spec[mask] = np.polyval(c[::-1], w)

    c = np.asarray([-1.025961e15,1.3172859e12,-3.873465e8,46486.541,-2.049])
    mask = (wavelength >= 3750) & (wavelength < 6250)
    w = wavelength[mask]    
    spec[mask] = np.polyval(c[::-1], w)

    c = np.asarray([4.861821e15,-2.2589885e12,4.3764376e8,-39279.61444,1.34388])
    mask = (wavelength >= 6250) & (wavelength < 8300)
    w = wavelength[mask]    
    spec[mask] = np.polyval(c[::-1], w)

    c = np.asarray([1.758394e15,-3.293986e11,1.6782617e7])
    mask = (wavelength >= 8300) & (wavelength < 8850)
    w = wavelength[mask]    
    spec[mask] = np.polyval(c[::-1], w)

    c = np.asarray([1.61455557e16,-6.544209e12,1.0159316e9,-70695.58136,1.852022])
    mask = (wavelength >= 8850) & (wavelength < 10000)
    w = wavelength[mask]    
    spec[mask] = np.polyval(c[::-1], w)

    c = np.asarray([7.97805136e14,-1.16906597e11,5.315222e6,-4.57327954,-3.473452e-3])
    mask = (wavelength >= 10000)
    w = wavelength[mask]    
    spec[mask] = np.polyval(c[::-1], w)

    if (units == 'frequency'):
        return spec * (wavelength * 1e-8)**2 / 299792458.0e2
    else:
        return spec

    return spec * (wavelength)

if (__name__ == '__main__'):
    import matplotlib.pyplot as pl
    x = np.linspace(3500, 8500, 200)
    res = contHSRA(x)
    pl.plot(x, res)