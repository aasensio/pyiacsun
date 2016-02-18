from __future__ import print_function

import numpy as np
import numpy.testing as npt
import pyiacsun
import time
import os

def test_voigt():
	voigt = pyiacsun.radtran.voigt(1.0, 0.1)
	npt.assert_equal(voigt, 0.37317014831126744)

def todo_test_milne():
	lambda0 = 6301.5080
	JUp = 2.0
	JLow = 2.0
	gUp = 1.5
	gLow = 1.833
	lambdaStart = 6300.8
	lambdaStep = 0.03
	nLambda = 100

	lineInfo = np.asarray([lambda0, JUp, JLow, gUp, gLow, lambdaStart, lambdaStep])

	wav = pyiacsun.radtran.milne.init(nLambda, lineInfo)

	lambda0 = 6302.5080
	JUp = 2.0
	JLow = 2.0
	gUp = 1.5
	gLow = 1.833
	lambdaStart = 6300.8
	lambdaStep = 0.03
	nLambda = 100

	lineInfo = np.asarray([lambda0, JUp, JLow, gUp, gLow, lambdaStart, lambdaStep])
	pyiacsun.radtran.milne.addLine(lineInfo)

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
		stokes = pyiacsun.radtran.milne.synth(nLambda, modelSingle, mu)
	end = time.time()
	print("Computing {0} models without derivatives took {1} s -> {2} models/s".format(nModels,(end-start),nModels / (end-start)))

def test_lte():
	data_path = os.path.join(pyiacsun.__path__[0], 'data')

	atmos = np.loadtxt(data_path+'/hsra_64.model', skiprows=2)
	lines = np.loadtxt(data_path+'/lines.dat')

	wl = np.linspace(10824,10830,200)
	pyiacsun.radtran.lte.initAtmos(atmos)
	pyiacsun.radtran.lte.initLines(lines, wl)

	start = time.time()
	for i in range(64):
		stokes, continuum = pyiacsun.radtran.lte.synthLines(atmos)
	print(time.time() - start)