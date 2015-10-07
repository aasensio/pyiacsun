milne
=====

Milne-Eddington synthesis

This is a module for the synthesis of spectral lines in the
Milne-Eddington approximation.

### Compilation
The code is written in Fortran 90 and has to be compiled with the gfortran free compiler.
Install it in your system and type the following in the `source` directory:

```
python setup.py build_ext --inplace
```

### Usage

```python
from milne import milne

# Define the details of the spectral line
lambda0 = 6301.5080
JUp = 2.0
JLow = 2.0
gUp = 1.5
gLow = 1.833
lambdaStart = 6300.8
lambdaStep = 0.03
nLambda = 50

lineInfo = np.asarray([lambda0, JUp, JLow, gUp, gLow, lambdaStart, lambdaStep])

# Define the model
BField = 100.0
BTheta = 20.0
BChi = 20.0
VMac = 2.0
damping = 0.0
B0 = 0.8
B1 = 0.2
mu = 1.0
VDop = 0.085
kl = 5.0
modelSingle = np.asarray([BField, BTheta, BChi, VMac, damping, B0, B1, VDop, kl])

# Generate the object and do one synthesis. The wavelength axis is a property of the object
s = milne(nLambda, lineInfo)
stokes = s.synth(modelSingle,mu)

# It can also synthesize many models at the same time
model = np.tile(modelSingle, (nModels,1)).T
stokes = s.synthGroup(model,mu)

# It can also return the response functions
stokes, dStokes = s.synthDerivatives(modelSingle,mu)

# Or for many models
stokes, dStokes = s.synthGroupDerivatives(model,mu)

```