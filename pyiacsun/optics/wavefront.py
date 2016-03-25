from math import *
import numpy as np
import sys
from scipy.special import gamma
from ipdb import set_trace as stop


"""
Routines for the computation of Zernike polynomials and seeing PSFs
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


def zernike(j,npix=256,phase=0.0):
  """
  Return the Zernike polynomial of Noll order j
  
  Args:
      j (int): Zernike order (Noll ordering)
      npix (int, optional): number of pixels of the final matrix
      phase (float, optional): phase
  
  Returns:
      real: circular phase with the output Zernike polynomial
  """
  if (j > 820):
    print("For n < 40, pick j < 820")
    sys.exit()

  x = np.arange(-npix/2,npix/2,dtype='d')
  y = np.arange(-npix/2,npix/2,dtype='d')

  xarr = np.outer(np.ones(npix,dtype='d'),x)
  yarr = np.outer(y,np.ones(npix,dtype='d'))

  rarr = np.sqrt(np.power(xarr,2) + np.power(yarr,2))/(npix/2)
  thetarr = np.arctan2(yarr,xarr) + phase

  outside = np.where(rarr > 1.0)

  narr = np.arange(40)
  jmax = (narr+1)*(narr+2)/2
  wh = np.where(j <= jmax)
  n = wh[0][0]
  mprime = j - n*(n+1)/2
  if ((n % 2) == 0):
    m = 2*int(floor(mprime/2))
  else:
    m = 1 + 2*int(floor((mprime-1)/2))

  radial = np.zeros((npix,npix),dtype='d')
  zarr = np.zeros((npix,npix),dtype='d')

  for s in range(int((n-m)/2 + 1)):
    tmp = pow(-1,s) * factorial(n-s)
    tmp /= factorial(s)*factorial((n+m)/2 - s)*factorial((n-m)/2 - s)
    radial += tmp*np.power(rarr,n-2*s)

  if (m == 0):
    zarr = radial
  else:
    if ((j % 2) == 0):
      zarr = sqrt(2.0)*radial*np.cos(m*thetarr)
    else:
      zarr = sqrt(2.0)*radial*np.sin(m*thetarr)
      m *= -1

  zarr *= sqrt(n+1)
  zarr[outside] = 0.0

  return zarr, n, m

def noll_to_nm(j):
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


def nollIndices(j):
  narr = np.arange(40)
  jmax = (narr+1)*(narr+2)/2
  wh = np.where(j <= jmax)
  n = wh[0][0]
  mprime = j - n*(n+1)/2
  if ((n % 2) == 0):
    m = 2*int(floor(mprime/2))
  else:
    m = 1 + 2*int(floor((mprime-1)/2))

  if ((j % 2) != 0):
    m *= -1

  return n, m

def aperture(npix=256, cent_obs=0.0, spider=0):
  """
  Compute the aperture image of a telescope
  
  Args:
      npix (int, optional): number of pixels of the aperture image
      cent_obs (float, optional): central obscuration fraction
      spider (int, optional): spider size in pixels
  
  Returns:
      real: returns the aperture of the telescope
  """
  illum = np.ones((npix,npix),dtype='d')
  x = np.arange(-npix/2,npix/2,dtype='d')
  y = np.arange(-npix/2,npix/2,dtype='d')

  xarr = np.outer(np.ones(npix,dtype='d'),x)
  yarr = np.outer(y,np.ones(npix,dtype='d'))

  rarr = np.sqrt(np.power(xarr,2) + np.power(yarr,2))/(npix/2)
  outside = np.where(rarr > 1.0)
  inside = np.where(rarr < cent_obs)

  illum[outside] = 0.0
  if np.any(inside[0]):
    illum[inside] = 0.0

  if (spider > 0):
   start = npix/2 - int(spider)/2
   illum[start:start+int(spider),:] = 0.0
   illum[:,start:start+int(spider)] = 0.0

  return illum

def plane_wave(npix=256):
  """
  Returns a plane wave image
  
  Args:
      npix (int, optional): number of pixels
  
  Returns:
      real: a plane wave
  """
  wf = np.zeros((npix,npix),dtype='d')

  return wf

def seeing(d_over_r0, npix=256, nterms=15, level=None, quiet=False):
  """
  Returns the wavefront resulting from a realization of the seeing
  
  Args:
      d_over_r0 (real): D/r0
      npix (int, optional): number of pixels
      nterms (int, optional): number of terms to include
      level (real, optional): precision level on the wavefront. This sets the number of terms
      quiet (bool, optional): verbose
  
  Returns:
      real: wavefront
  """
  scale = pow(d_over_r0,5.0/3.0)

  if level:
    narr = np.arange(400,dtype='d') + 2
    coef = np.sqrt(0.2944*scale*(np.power((narr-1),-0.866) - np.power(narr,-0.866)))
    wh = np.where(coef < level)
    n = wh[0][0]
    norder = int(ceil(sqrt(2*n)-0.5))
    nterms = norder*(norder+1)/2
    if (nterms < 15):
      nterms = 15

  wf = np.zeros((npix,npix),dtype='d')

  if (nterms == 0):
    return wf

  resid = np.zeros(nterms,dtype='d')
  coeff = np.zeros(nterms,dtype='d')

  resid[0:10] = [1.030,0.582,0.134,0.111,0.088,0.065,0.059,0.053,0.046,0.040]
  if (nterms > 10):
    for i in range(10,nterms):
      resid[i] = 0.2944*pow(i+1,-0.866)

  for j in range(2,nterms+1):
    coeff[j-1] = sqrt((resid[j-2]-resid[j-1])*scale)
    wf += coeff[j-1]*np.random.normal()*zernike(j,npix=npix)

  if not quiet:
    print("Computed Zernikes to term %d and RMS %f" % (nterms,coeff[nterms-1]))

  return wf

def psf(aperture, wavefront, overfill=1):
  """
  Transform an aperture and the wavefront to a PSF
  
  Args:
      aperture (TYPE): aperture
      wavefront (TYPE): wavefront
      overfill (int, optional): number of extra pixels
  
  Returns:
      real: PSF
  """
  npix = len(wavefront)
  nbig = npix*overfill
  wfbig = np.zeros((nbig,nbig),dtype='d')

  half = (nbig - npix)/2
  wfbig[half:half+npix,half:half+npix] = wavefront

  illum = np.zeros((nbig,nbig),dtype='d')
  illum[half:half+npix,half:half+npix] = aperture

  phase = np.exp(wfbig*(0.+1.j))
  input = illum*phase

  ft = np.fft.fft2(input)
  powft = np.real(np.conj(ft)*ft)

  sorted = np.zeros((nbig,nbig),dtype='d')
  sorted[:nbig/2,:nbig/2] = powft[nbig/2:,nbig/2:]
  sorted[:nbig/2,nbig/2:] = powft[nbig/2:,:nbig/2]
  sorted[nbig/2:,:nbig/2] = powft[:nbig/2,nbig/2:]
  sorted[nbig/2:,nbig/2:] = powft[:nbig/2,:nbig/2]

  crop =  sorted[half:half+npix,half:half+npix]

  # fluxrat = np.sum(crop)/np.sum(sorted)
  #print "Cropped PSF has %.2f%% of the flux" % (100*fluxrat)

  return crop

def psfScale(D, wavelength, pixSize):
  """
  Return the PSF scale appropriate for the required pixel size, wavelength and telescope diameter
  The aperture is padded by this amount; resultant pix scale is lambda/D/psf_scale, so for instance full frame 256 pix
  for 3.5 m at 532 nm is 256*5.32e-7/3.5/3 = 2.67 arcsec for psf_scale = 3
  
  Args:
      D (real): telescope diameter in m
      wavelength (real): wavelength in Angstrom
      pixSize (real): pixel size in arcsec
  
  Returns:
      real: psf scale
  """
  DInCm = D * 100.0
  wavelengthInCm = wavelength * 1e-8
  return 206265.0 * wavelengthInCm / (DInCm * pixSize)

if __name__ == "__main__":
  Z, n, m = zernike(1)
  gradZ = np.gradient(Z)

  # N = 256           # number of pixels across
  # D_over_r0 = 20.0        # seeing is  lambda/r0 radians

  # illum = aperture(npix=N, cent_obs = 0.3,spider=2)
  # wf = seeing(D_over_r0, npix=N, nterms=20)

  # psf_scale = 3 # pads aperture by this amount; resultant pix scale is
  #     # lambda/D/psf_scale, so for instance full frame 256 pix
  #     # for 3.5 m at 532 nm is 256*5.32e-7/3.5/3 = 2.67 arcsec
  #     # for psf_scale = 3

  # D = 1.5
  # wavelength = 6301.0
  # pixSize = 23.71 / 725.0
  # # generate speckle pattern given my wavefront and aperture illumination
  # psfFinal = psf(illum,wf,overfill=psfScale(D, wavelength, pixSize))