import pymilne_i686 as pymilne
import numpy as np
import matplotlib.pyplot as pl
import time

lambda0 = 6301.5080
JUp = 2.0
JLow = 2.0
gUp = 1.5
gLow = 1.833
lambdaStart = 6300.8
lambdaStep = 0.03
nLambda = 100

lineInfo = np.asarray([lambda0, JUp, JLow, gUp, gLow, lambdaStart, lambdaStep])

wav = pymilne.init(nLambda, lineInfo)

lambda0 = 6302.5080
JUp = 2.0
JLow = 2.0
gUp = 1.5
gLow = 1.833
lambdaStart = 6300.8
lambdaStep = 0.03
nLambda = 100

lineInfo = np.asarray([lambda0, JUp, JLow, gUp, gLow, lambdaStart, lambdaStep])
pymilne.addLine(lineInfo)

stokes = np.zeros((4,nLambda))

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

nModels = 5000
start = time.time()
for i in range(nModels):
	stokes = pymilne.synth(nLambda, modelSingle, mu)
end = time.time()
print "Computing {0} models without derivatives took {1} s -> {2} models/s".format(nModels,(end-start),nModels / (end-start))

# model = np.tile(modelSingle, (nModels,1)).T
# start = time.time()
# stokes = s.synthGroup(model,mu)
# end = time.time()
# print "Computing {0} models without derivatives took {1} s -> {2} models/s".format(nModels,(end-start),nModels / (end-start))

# start = time.time()
# for i in range(nModels):
# 	stokes, dStokes = s.synthDerivatives(modelSingle,mu)
# end = time.time()
# print "Computing {0} models with derivatives took {1} s -> {2} models/s".format(nModels,(end-start),nModels / (end-start))

# model = np.tile(modelSingle, (nModels,1)).T
# start = time.time()
# stokes = s.synthGroupDerivatives(model,mu)
# end = time.time()
# print "Computing {0} models with derivatives took {1} s -> {2} models/s".format(nModels,(end-start),nModels / (end-start))