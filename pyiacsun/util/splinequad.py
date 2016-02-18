__all__ = ['splinequad']

import numpy as np

def splinequad(x):
    """Return the quadrature weights when a function is interpolated with cubic splines
    The function retuns the values of the weights for a function that is evaluated in the locations x,
    so that

    int(f(x),x=a..b) = sum(w*f(x)), with x a partition of the interval [a,b]

    It is based on the formulae developed in "The use of quadrature weights in cubic spline integration", C. M. Leung, R. W. Quan
    International Journal of Mathematical Education in Science and Technology, Volume 15, Issue 3, 1984
    http://www.tandfonline.com/doi/abs/10.1080/0020739840150306
    
    Args:
        x (array): abscissas where the function is evaluated
    
    Returns:
        float: array with the weights
    """
    n = len(x)
    dx = np.zeros(n-1)
    for i in range(n-1):
        dx[i] = x[i+1] - x[i]

    w1 = np.zeros(n)
    w2 = np.zeros(n)
    a = np.zeros(n)
    b = np.zeros(n)
    c = np.zeros(n)
    alpha = np.zeros(n)
    beta = np.zeros(n)
    gamma = np.zeros(n)
    
    w1[0] = 0.5 * dx[0]
    w1[n-1] = 0.5 * dx[n-2]

    w2[0] = -1.0 / 24.0 * dx[0]**3
    w2[n-1] = -1.0 / 24.0 * dx[n-2]**3

    a[0] = -1.0 / dx[0]
    a[n-1] = -1.0 / dx[n-3]

    c[0] = -1.0 / dx[1]
    c[n-1] = -1.0 / dx[n-2]

    b[0] = -(a[0]+c[0])
    a[n-1] = -(a[n-1] + c[n-1])
    for i in range(1,n-1):
        w1[i] = 0.5 * ( dx[i-1] + dx[i] )

        w2[i] = -1.0 / 24.0 * ( dx[i-1]**3 + dx[i]**3 )

        a[i] = dx[i-1]
        c[i] = dx[i]
        b[i] = 2.0*(a[i] + c[i])

        alpha[i] = 6.0 / dx[i-1]
        gamma[i] = 6.0 / dx[i]
        beta[i] = -(alpha[i] + gamma[i])


    U = np.zeros((n,n))
    T = np.zeros((n,n))

    U[0,0:3] = np.asarray([a[0],b[0],c[0]])
    U[n-1,-3:] = np.asarray([a[n-1],b[n-1],c[n-1]])
   
    for i in range(1,n-1):
        U[i,i-1:i+2] = np.asarray([a[i], b[i], c[i]])
        T[i,i-1:i+2] = np.asarray([alpha[i], beta[i], gamma[i]])

    return w1 + w2.dot(np.linalg.inv(U).dot(T))