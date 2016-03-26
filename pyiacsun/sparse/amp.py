import scipy.linalg
import numpy as np
import matplotlib.pyplot as pl
from ipdb import set_trace as stop

def eta(x, threshold):
    return np.sign(x) * np.fmax(np.abs(x) - threshold, 0)

def etaprime(x, threshold):
    return (x > threshold) + (x < -threshold)

def largestElement(x, n):
    lenx = len(x)
    if (n > lenx):
        n = lenx-1
    if (n < 0):
        n = 0

    t = np.sort(x)[::-1]
    return t[n]

def amp(A, AT, x0, eta, etaprime, b, mu, maxIter, xOK, tol):
    xhat = np.copy(x0)
    z = np.copy(b)
    
    delta = 1.0 * len(b) / len(xhat)

    err = np.zeros(maxIter)

    gamm = 0.0

    for i in range(maxIter):
        xhat = eta(xhat + AT(z), mu + gamm)

        z = b - A(xhat) + z / delta * np.mean(etaprime(xhat + AT(z), mu + gamm))

        gamm = (mu+gamm) / delta * np.mean(etaprime(xhat + AT(z), mu + gamm))

        stopping = np.linalg.norm(b - A(xhat)) / np.linalg.norm(b)
        err[i] = stopping
        print (i, stopping, mu+gamm)

        if (stopping < tol):
            return xhat, err
    return xhat, err

if (__name__ == "__main__"):
    M = 200
    N = 1000
    K = 10

    mu = 0.005
    sigma = 0.00001

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

    sol, err = amp(A, At, x0, eta, etaprime, b, mu, 500, x, 1e-6)

    pl.plot(sol)
    pl.plot(x, 'o')