lte
===

Routines for the synthesis of Stokes profiles in LTE. The routines are
written in IDL, Fortran 90 and Python.

## python

Python wrapper for the synthesis of intensity profiles neglecting polarization. The code needs
to be first compiler using an appropriate Fortran compiler. I have tested it with GFortran.
To compile the code, just enter the `source` directory and type

	python setup.py build_ext --inplace

This will generate the pylte.so library that is later imported from Python. The code is run using
the following calls
```python
import pylte
pylte.initAtmos(atmos)
pylte.initLines(lines, wl)
```

where `atmos` is an array of shape [ndepth,4] that contains the following columns

	- log(tau) 5000 Angstrom
	- Temperature [K]
	- microturbulence velocity [cm/s]
	- macroscopic velocity [km/s]

and `lines` is an array of shape [nlines,11] containing the following information that
describes the spectral lines to be synthesized:

	- Central wavelength [Angstrom]
	- Atomic number of the species
	- Ionization state (1 is neutral, 2 is one time ionized, etc.)
	- log(gf)
	- Energy of the lower level in cm^-1
	- Landé factor of the upper level
	- Landé factor of the lower level
	- Total angular momentum J of the upper level
	- Total angular momentum J of the upper level
	- sigma value of the Barklem & O'Mara broadening (set to zero to avoid using it)
	- alpha value of the Barklem & O'Mara broadening (set to zero to avoid using it)

Finally, `wl` is the wavelength axis defined as a vector.

The emergent Stokes I profile can be obtained by just calling

```python
out = pylte.synthLines(atmos)
```

The last column of `out` contains the continuum.

## pythonZeeman

Python wrapper for the synthesis of intensity profiles neglecting polarization. The code needs
to be first compiler using an appropriate Fortran compiler. I have tested it with GFortran.
To compile the code, just enter the `source` directory and type

	python setup.py build_ext --inplace

This will generate the pylte.so library that is later imported from Python. The code is run using
the following calls
```python
import pylte
pylte.initAtmos(atmos)
pylte.initLines(lines, wl)
```

where `atmos` is an array of shape [ndepth,4] that contains the following columns

	- log(tau) 5000 Angstrom
	- Temperature [K]
	- microturbulence velocity [cm/s]
	- macroscopic velocity [km/s]
	- Magnetic field strength [G]
	- Magnetic field inclination [deg]
	- Magnetic field azimuth [deg]

and `lines` is an array of shape [nlines,11] containing the following information that
describes the spectral lines to be synthesized:

	- Central wavelength [Angstrom]
	- Atomic number of the species
	- Ionization state (1 is neutral, 2 is one time ionized, etc.)
	- log(gf)
	- Energy of the lower level in cm^-1
	- Landé factor of the upper level
	- Landé factor of the lower level
	- Total angular momentum J of the upper level
	- Total angular momentum J of the upper level
	- sigma value of the Barklem & O'Mara broadening (set to zero to avoid using it)
	- alpha value of the Barklem & O'Mara broadening (set to zero to avoid using it)

Finally, `wl` is the wavelength axis defined as a vector.

The emergent Stokes I profile can be obtained by just calling

```python
stokes, continuum = pylte.synthLines(atmos)
```