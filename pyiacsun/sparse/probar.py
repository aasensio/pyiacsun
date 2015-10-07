import pyFasta
import scipy.linalg
import numpy as np
import matplotlib.pyplot as pl
from ipdb import set_trace as stop

M = 200
N = 1000
K = 10
mu = 0.02
sigma = 0.01

# Create sparse signal
x = np.zeros(N)
ind = np.random.permutation(N)
x[ind[0:K]] = 1.0

# Define matrix
AMat = np.random.normal(size=(M,N))
AMat /= np.linalg.norm(AMat, 2)

# Define observation vector
b = AMat.dot(x)
b += np.random.normal(scale=sigma, size=b.shape)

# Initial state
x0 = np.zeros(N)

A = lambda z : AMat.dot(z)
At = lambda z : AMat.T.dot(z)

f = lambda z: 0.5*np.linalg.norm(z-b)**2
fgrad = lambda z : z-b

h = lambda z : mu * np.linalg.norm(z,1)
hprox = lambda z, t : pyFasta.proxes.prox_l1(z, t*mu)

out = pyFasta.fasta(A, At, f, fgrad, h, hprox, x0, backtrack=True, verbose=True)
res = out.optimize()

pl.stem(x, 'r')
pl.stem(res, 'b')