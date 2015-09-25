from milne import milne as milne
import numpy as np
import matplotlib.pyplot as pl
import datetime as dt

lambda0 = 23000.0
JUp = 3.0
JLow = 2.0
gUp = 1.5 / 1800
gLow = 1.5 / 1800
lambdaStart = 22999
lambdaStep = 0.01
nLambda = 200

lineInfo = [lambda0, JUp, JLow, gUp, gLow, lambdaStart, lambdaStep, nLambda]

s = milne(lineInfo)

stokes = np.zeros((4,nLambda))

BField = [0.0,500.0,1000.0,1500.0,5000.,10000.,1e5,1e6]
BTheta = 20.0
BChi = 20.0
VMac = 0.0
damping = 0.0
beta = 3.0
mu = 1.0
VDop = 0.15
kl = 5.0

f = pl.figure()
ax = f.add_subplot(1,1,1)
for i in range(len(BField)):
	model = [BField[i], BTheta, BChi, VMac, damping, beta, VDop, kl]
	wavelength, stokes = s.synth(model,mu)
	ax.plot(wavelength-lambda0, stokes[0,:], label=str(BField[i])+' G')

pl.legend(loc='lower right')	
f.tight_layout()
f.savefig("HF.pdf")
