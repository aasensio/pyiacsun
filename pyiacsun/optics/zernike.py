import numpy as np
from scipy.special import gamma
#from ipdb import set_trace as stop
#import matplotlib.pyplot as pl

__all__ = ['zernikeMachine']

"""
Routines for the computation of Zernike polynomials 
http://www-physics.ucsd.edu/~tmurphy/astr597/exercises/speckle
"""

def factorial(n):
    """
    Compute the factorial of a number
  
    Args:
      n (real): input number
  
    Returns:
      real: n!
    """
    return gamma(n+1)

def iseven(n):
    return n % 2 == 0

class zernikeMachine(object):
    """
    Class that can be used to carry out several operations with Zernike functions
    """
    def __init__(self, npix=256, phase=0.0):
        """Class that can be used to carry out several operations with Zernike functions
        
        Args:
            npix (int, optional): number of pixels in which the Zernike functions will be returned
            phase (float, optional): additional phase to be added
        """
        self.npix = npix
        x = np.arange(-self.npix/2,self.npix/2,dtype='d')
        y = np.arange(-self.npix/2,self.npix/2,dtype='d')

        xarr = np.outer(np.ones(self.npix,dtype='d'),x)
        yarr = np.outer(y,np.ones(self.npix,dtype='d'))

        self.rarr = np.sqrt(np.power(xarr,2) + np.power(yarr,2))/(npix/2)
        self.thetarr = np.arctan2(yarr,xarr) + phase

        self.outside = np.where(self.rarr > 1.0)

        self.radial = np.zeros((self.npix,self.npix),dtype='d')
        self.zarr = np.zeros((self.npix,self.npix),dtype='d')

    def noll_to_nm(self, j):
        """Transform Noll indices to n,m pairs Z_j -> Z_n^m
        
        Args:
            j (int): Noll index
        
        Returns:
            int: (n,m) indices of the Zernike functions
        """
        n = int(np.floor(np.sqrt(2.0*j-0.75)-0.50)  )
        nn=(n-1)*(n+2)/2
        m=j-nn-2
        if np.ceil(n/2.) != np.floor(n/2.):
          m = 2*int(m/2.)+1
        else:
          m=2*int((m+1)/2.)
        if np.ceil(j/2.) != np.floor(j/2.):
          m=-m
        return n,m

    def nm_to_noll(n, m):
        """Transform n,m indices to Noll indices Z_n^m -> Z_j
        
        Args:
            n (int): n-index
            m (Tint): m-index
        
        Returns:
            int: Noll index of the Zernike functions
        """
        signo = 1
        if (m != 0):
            signo = m / np.abs(m)
        noll = 0

        if (np.abs(m) <= n):
            num_m = np.abs(m)
            j = (n-1)*(n+2)/2.0+2
            if (np.floor((n-m)/2.) == np.ceil((n-m)/2.)):
                noll = j + num_m
            if (m != 0):
                if (signo == -1 and np.floor(noll/2.0) == np.ceil(noll/2.0)):
                    noll -= 1
                if (signo == 1 and np.floor(noll/2.0) != np.ceil(noll/2.0)):
                    noll -= 1
        return int(noll)

    def zernikeNoll(self, j):
        """
        Return the Zernike polynomial of Noll order j
      
        Args:
            j (int): Zernike order (Noll ordering)
              
        Returns:
            real: circular phase with the output Zernike polynomial
        """
        if (j > 820):
            print("For n < 40, pick j < 820")
            return 0

        n, m = self.noll_to_nm(j)

        return self.zernikeNM(n, m)

    def zernikeNM(self, n, m):
        """
        Return the Zernike polynomial of indices n and m
      
        Args:
            n, m (int): orders
              
        Returns:
            real: circular phase with the output Zernike polynomial
        """

        self.radial *= 0

        for s in range(int((n-np.abs(m))/2 + 1)):
            tmp = pow(-1,s) * factorial(n-s)
            tmp /= factorial(s)*factorial((n+np.abs(m))/2.0 - s)*factorial((n-np.abs(m))/2.0 - s)
            self.radial += tmp*np.power(self.rarr,n-2*s)

        if (m == 0):
            self.zarr = np.copy(self.radial)
        else:
            if (m >= 0):
                self.zarr = np.sqrt(2.0)*self.radial*np.cos(np.abs(m)*self.thetarr)
            else:
                self.zarr = np.sqrt(2.0)*self.radial*np.sin(np.abs(m)*self.thetarr)

        self.zarr *= np.sqrt(n+1)

        self.zarr[self.outside] = 0.0

        return self.zarr

    def gradZernike(self, j):
        """X and Y gradients of the Zernike function of index j
        
        Args:
            j (int): Noll index
        
        Returns:
            float: dZ/dX and dZ/dY
        """
        n, m = self.noll_to_nm(j)
        alpham = 1.0
        if (m < 0):
            alpham = -1.0

        dZernX = np.zeros_like(self.thetarr)
        dZernY = np.zeros_like(self.thetarr)
        for nprime in range(n):
            temp1 = self._aFactor(n, m, nprime, m-1) * self.zernikeNM(nprime, alpham * np.abs(m-1))
            temp2 = np.sign(m+1) * alpham * self._aFactor(n, m, nprime, m+1) * self.zernikeNM(nprime, alpham * np.abs(m+1))
            dZernX += temp1 + temp2

            temp1 = -np.sign(m-1) * alpham * self._aFactor(n, m, nprime, m-1) * self.zernikeNM(nprime, -alpham * np.abs(m-1))
            temp2 = self._aFactor(n, m, nprime, m+1) * self.zernikeNM(nprime, -alpham * np.abs(m+1))
            dZernY += temp1 + temp2

        return dZernX, dZernY
        
    def _aFactor(self, n, m, nprime, mp):
        out = 0.0
        if (n >= np.abs(m) and iseven(n-np.abs(m)) and nprime >= np.abs(mp) and iseven(nprime - np.abs(mp))):
            num = 2.0
            den = 2.0
            if (m == 0):
                num = 1.0
            if (mp == 0):
                den = 1.0
            out = np.sqrt(num / den * (n+1) * (nprime+1))
        return out

    def generateTurbulentZernikesKolmogorov(self, D_over_r0, nModes, firstMode=1):
        """Generate the covariance matrix for the Zernike coefficients for a given value of r0 using Kolmogorov statistics
        
        Args:
            r0 (float): Fried radius [m]
        
        Returns:
            None
        """
        covariance = np.zeros((nModes,nModes))
        for i in range(nModes):
            ni, mi = self.noll_to_nm(i+firstMode)
            for j in range(nModes):
                nj, mj = self.noll_to_nm(j+firstMode)
                if (iseven(i - j)):
                    if (mi == mj):
                        phase = (-1.0)**(0.5*(ni+nj-2*mi))
                        t1 = np.sqrt((ni+1)*(nj+1)) * np.pi**(8.0/3.0) * 0.0072 * (D_over_r0)**(5.0/3.0)
                        t2 = gamma(14./3.0) * gamma(0.5*(ni+nj-5.0/3.0))
                        t3 = gamma(0.5*(ni-nj+17.0/3.0)) * gamma(0.5*(nj-ni+17.0/3.0)) * gamma(0.5*(ni+nj+23.0/3.0))
                        covariance[i,j] = phase * t1 * t2 / t3

        return np.random.multivariate_normal(np.zeros(nModes), covariance)


if (__name__ == "__main__"):
    import matplotlib.pyplot as pl
    zern = zernikeMachine()
    Z = zern.zernikeNoll(4)
    dX, dY = zern.gradZernike(4)
    f, ax = pl.subplots(ncols=3, nrows=1, figsize=(15,6))
    ax[0].imshow(Z)
    ax[1].imshow(dX)
    ax[2].imshow(dY)
