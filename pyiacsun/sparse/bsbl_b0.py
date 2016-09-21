from __future__ import print_function

__all__ = ["bsbl_b0"]

import numpy as np
import copy
import scipy.special as sp
import scipy.linalg as sl
import matplotlib.pyplot as pl

def bsbl_b0(Phi, Y, blkStartLoc, LearnLambda, learntype=1, prune_gamma=0.001, 
    lambdaPar=1e-14, max_iters=600, epsilon=1e-8, verbose=False):
    """
    Recover block sparse signal (1D) exploiting intra-block correlation, given the block partition.

          The algorithm solves the inverse problem for the block sparse
          model with known block partition:
                        y = Phi * x + v
    
    Args:
        Phi (float): N X M known matrix
        Y (float): N X 1 measurement vector 
        blkStartLoc (int): Start location of each block
        LearnLambda (int): (1) If LearnLambda = 1, use the lambda learning rule for very LOW SNR cases (SNR<10dB)
                               (using lambda=std(y)*1e-2 or user-input value as initial value)
                           (2) If LearnLambda = 2, use the lambda learning rule for medium noisy cases (SNR>10dB) 
                               (using lambda=std(y)*1e-2 or user-input value as initial value)
                           (3) If LearnLambda = 0, do not use the lambda learning rule 
                               ((using lambda=1e-14 or user-input value as initial value)
        learntype (int, optional): LEARNTYPE = 0: Ignore intra-block correlation
                                   LEARNTYPE = 1: Exploit intra-block correlation 
        prune_gamma (float, optional): threshold to prune out small gamma_i 
                                       (generally, 10^{-3} or 10^{-2})
        lambdaPar (float, optional): user-input value for lambda
                                    [ Default: LAMBDA=1e-14 when LearnLambda=0; LAMBDA=std(y)*1e-2 in noisy cases]
        max_iters (int, optional): Maximum number of iterations.
        epsilon (float, optional): Solution accuracy tolerance parameter 
        verbose (bool, optional): Display flag. If True: show output; If False: suppress output
    
    Returns:
        dict: Result  
            Result.x          : the estimated block sparse signal
            Result.gamma_used : indexes of nonzero groups in the sparse signal
            Result.gamma_est  : the gamma values of all the groups of the signal
            Result.B          : the final value of the B
            Result.count      : iteration times

    Examples:

    For most noisy environment (SNR > 10dB):
          
           Result =  BSBL_BO(Phi, y, blkStartLoc, 2);  

    For very low SNR cases (SNR < 10 dB):
           
           Result =  BSBL_BO(Phi, y, blkStartLoc, 1);   

    For noiseless cases:
          
           Result =  BSBL_BO(Phi, y, blkStartLoc, 0); 

    Author:

    [1] Zhilin Zhang, Bhaskar D. Rao, Extension of SBL Algorithms for the 
       Recovery of Block Sparse Signals with Intra-Block Correlation, 
       available at: http://arxiv.org/abs/1201.0862
  
    [2] Zhilin Zhang, Tzyy-Ping Jung, Scott Makeig, Bhaskar D. Rao, 
       Low Energy Wireless Body-Area Networks for Fetal ECG Telemonitoring 
       via the Framework of Block Sparse Bayesian Learning, 
       available at: http://arxiv.org/pdf/1205.1287v1.pdf

    [3] webpage: http://dsp.ucsd.edu/~zhilin/BSBL.html, or
                 https://sites.google.com/site/researchbyzhang/bsbl
 
    Zhilin Zhang (z4zhang@ucsd.edu, zhangzlacademy@gmail.com)
    """

    scl = np.std(Y)
    if ((scl < 0.4) or (scl > 1)):
        Y = Y / scl * 0.4

    if (LearnLambda == 0):
        prune_gamma = 1e-3
    elif (LearnLambda == 2):
        lambdaPar = scl * 1e-2
        prune_gamma = 1e-2
    elif (LearnLambda == 1):
        lambdaPar = scl * 1e-2
        prune_gamma = 1e-2
    else:
        print("Unrecognized value for LearnLambda")
        exit()

    if (verbose):
        print('====================================================')
        print('           Running BSBL-BO ....... ')
        print('           Information about parameters...')
        print('====================================================')
        print('PRUNE_GAMMA  : {0}'.format(prune_gamma))
        print('lambda       : {0}'.format(lambdaPar))
        print('LearnLambda  : {0}'.format(LearnLambda))
        print('LearnType    : {0}'.format(learntype))
        print('EPSILON      : {0}'.format(epsilon))
        print('MAX_ITERS    : {0}'.format(max_iters))

# Initialization
    N, M = Phi.shape
    Phi0 = np.copy(Phi)
    blkStartLoc0 = copy.deepcopy(blkStartLoc)
    p = len(blkStartLoc)
    blkLenList = []
    
    for k in range(p-1):
        blkLenList.append(blkStartLoc[k+1] - blkStartLoc[k])

    blkLenList.append(M - blkStartLoc[-1])
    blkLenList = np.array(blkLenList)

    maxLen = np.max(blkLenList)
    if (np.sum(blkLenList == maxLen) == p):
        equalSize = 1
    else:
        equalSize = 0

    Sigma0 = [None] * p
    for k in range(p):
        Sigma0[k] = np.eye(blkLenList[k])

    gamma = np.ones((p,1))
    keep_list = np.arange(p).T
    usedNum = len(keep_list)
    mu_x = np.zeros((M,1))
    count = 0

    while(True):
        count += 1

# =========== Prune weights as their hyperparameters go to zero ==============
        if (np.min(gamma) < prune_gamma):
            index = np.where(gamma > prune_gamma)[0]
            usedNum = len(index)
            keep_list = keep_list[index] 
            if (keep_list.size == 0):
                print('\n====================================================================================\n');
                print('x becomes zero vector. The solution may be incorrect. \n')
                print('Current ''prune_gamma'' = {0}, and Current ''EPSILON'' = {1}.\n'.format(prune_gamma,epsilon))
                print('Try smaller values of ''prune_gamma'' and ''EPSILON'' or normalize ''y'' to unit norm.\n')
                print('====================================================================================\n\n')
                break

            blkStartLoc = blkStartLoc[index]
            blkLenList = blkLenList[index]

# Prune gamma and associated components in Sigma0
            gamma = gamma[index]
            temp = copy.deepcopy(Sigma0)
            Sigma0 = [None] * usedNum
            for k in range(usedNum):
                Sigma0[k] = temp[index[k]]

# Construct new Phi
            temp = []
            for k in range(usedNum):                
                temp.append(Phi0[:,blkStartLoc[k]:blkStartLoc[k]+blkLenList[k]])

            Phi = np.hstack(temp)

# =========== Compute new weights ==============
        mu_old = np.copy(mu_x)
        PhiBPhi = np.zeros((N,N))
        currentLoc = -1
        for i in range(usedNum):
            currentLen = Sigma0[i].shape[0]
            currentLoc += 1
            currentSeg = np.arange(currentLoc, currentLoc + currentLen)

            PhiBPhi += Phi[:,currentSeg] @ Sigma0[i] @ Phi[:,currentSeg].T

            currentLoc = currentSeg[-1]

        num = Phi.T
        den = PhiBPhi + lambdaPar * np.eye(N)
        H = np.linalg.lstsq(den.T, num.T)[0].T
        Hy = H @ Y
        HPhi = H @ Phi

        mu_x = np.zeros((Phi.shape[1],1))
        Sigma_x = [None] * usedNum
        Cov_x = [None] * usedNum
        B = [None] * usedNum
        invB = [None] * usedNum
        B0 = np.zeros((maxLen,maxLen))
        r0 = np.zeros(1)
        r1 = np.zeros(1)
        currentLoc = -1

        for i in range(usedNum):
            currentLen = Sigma0[i].shape[0]
            currentLoc += 1
            seg = np.arange(currentLoc, currentLoc + currentLen)

            mu_x[seg] = Sigma0[i] @ Hy[seg]
            Sigma_x[i] = Sigma0[i] - Sigma0[i] @ HPhi[np.ix_(seg,seg)] @ Sigma0[i]
            Cov_x[i] = Sigma_x[i] + mu_x[seg] @ mu_x[seg].T

            currentLoc = seg[-1]

            #=========== Learn correlation structure in blocks ===========
            # do not consider correlation structure in each block
            if (learntype == 0):
                B[i] = np.eye(currentLen)
                invB[i] = np.eye(currentLen)

            # constrain all the blocks have the same correlation structure
            elif (learntype == 1):
                if (equalSize == 0):
                    if (currentLen > 1):
                        temp = Cov_x[i] / gamma[i]
                        r0 += np.mean(np.diag(temp))
                        r1 += np.mean(np.diag(temp))
                        # CHECK!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                elif (equalSize == 1):
                    B0 += Cov_x[i] / gamma[i]

        
#=========== Learn correlation structure in blocks with Constraint 1 ===========
# If blocks have the same size
        if ((equalSize == 1) and (learntype == 1)):

            # Constrain all the blocks have the same correlation structure
            # (an effective strategy to avoid overfitting)

            b = np.mean(np.diag(B0, k=1)) / np.mean(np.diag(B0))

            if (np.abs(b) >= 0.99):
                b = 0.99 * np.sign(b)
            bs = np.zeros(maxLen)
            for j in range(maxLen):
                bs[j] = b**j

            B0 = sl.toeplitz(bs)

            for i in range(usedNum):
                B[i] = B0
                invB[i] = np.linalg.inv(B[i])

        elif ((equalSize == 0) and (learntype == 1)):
            pass

            # FINISH!!!!!!!!!!!!!!!!!!!!!!!

        # estimate gamma(i) and lambda 
        if (LearnLambda == 1):
            gamma_old = np.copy(gamma)
            lambdaComp = 0
            currentLoc = -1
            for i in range(usedNum):
                currentLen = Sigma_x[i].shape[0]
                currentLoc += 1
                seg = np.arange(currentLoc, currentLoc + currentLen)

                gamma[i] = gamma_old[i] * np.linalg.norm(sl.sqrtm(B[i]) @ Hy[seg]) / np.sqrt(np.trace(HPhi[np.ix_(seg,seg)] @ B[i]))

                lambdaComp += np.trace( Phi[:,seg] @ Sigma_x[i] @ Phi[:,seg].T)

                Sigma0[i] = B[i] * gamma[i]

                currentLoc = seg[-1]

            lambdaPar = np.linalg.norm(Y - Phi @ mu_x, 2)**2 / N + lambdaComp / N

                # gamma[i] = gamma_old[i] * 
        elif (LearnLambda == 2):
            gamma_old = np.copy(gamma)
            lambdaComp = 0
            currentLoc = -1
            for i in range(usedNum):
                currentLen = Sigma_x[i].shape[0]
                currentLoc += 1
                seg = np.arange(currentLoc, currentLoc + currentLen)

                gamma[i] = gamma_old[i] * np.linalg.norm(sl.sqrtm(B[i]) @ Hy[seg]) / np.sqrt(np.trace(HPhi[np.ix_(seg,seg)] @ B[i]))

                lambdaComp += np.trace( Sigma_x[i] @ invB[i]) / gamma_old[i]

                Sigma0[i] = B[i] * gamma[i]

                currentLoc = seg[-1]

            lambdaPar = np.linalg.norm(Y - Phi @ mu_x, 2)**2 / N + lambdaPar * (len(mu_x) - lambdaComp) / N
        else:
            gamma_old = np.copy(gamma)
            lambdaComp = 0
            currentLoc = -1
            for i in range(usedNum):
                currentLen = Sigma_x[i].shape[0]
                currentLoc += 1
                seg = np.arange(currentLoc, currentLoc + currentLen)

                gamma[i] = gamma_old[i] * np.linalg.norm(sl.sqrtm(B[i]) @ Hy[seg]) / np.sqrt(np.trace(HPhi[np.ix_(seg,seg)] @ B[i]))

                Sigma0[i] = B[i] * gamma[i]

                currentLoc = seg[-1]

            
        # ================= Check stopping conditions, etc. ==============
        if (len(mu_x) == len(mu_old)):
            dmu = np.max(np.abs(mu_old - mu_x))
            if (dmu < epsilon):
                break
        if (verbose):
            print("Iters: {0:3d} - Num. coeffs: {1:3d} - minGamma: {2:12.8f} - gammaChange: {3:12.8f} - muChange: {4:12.8f}".format(count, usedNum, np.min(gamma), np.max(np.abs(gamma-gamma_old)), dmu))

        if (count >= max_iters):
            print("Reached max iterations. Stop")
            break

    if (keep_list.size == 0):
        result = {}
        result['x'] = np.zeros((M,1))
        result['gamma_used'] = []
        result['gamma_est'] = np.zeros((p,1))
        result['B'] = B
        result['count'] = count
        result['lambdatrace'] = lambdaPar

    else:
        result = {}

        # Expand hyperparameters
        gamma_used = np.sort(keep_list)
        gamma_est = np.zeros((p,1))            
        gamma_est[keep_list,:] = gamma

        # Reconstruct the original signal
        x = np.zeros((M,1))            
        currentLoc = -1
        for i in range(usedNum):
            currentLen = Sigma_x[i].shape[0]
            currentLoc += 1
            seg = np.arange(currentLoc, currentLoc + currentLen)

            realLocs = np.arange(blkStartLoc0[keep_list[i]], blkStartLoc0[keep_list[i]]+currentLen)

            x[realLocs] = mu_x[seg]

            currentLoc = seg[-1]

        if ((scl < 0.4) or (scl > 1)):
            result['x'] = x * scl / 0.4
        else:
            result['x'] = x
        result['gamma_used'] = gamma_used
        result['gamma_est'] = gamma_est
        result['B'] = B
        result['count'] = count
        result['lambdatrace'] = lambdaPar

        return result    
                


def randn2(*args,**kwargs):
    '''
    Calls rand and applies inverse transform sampling to the output.
    '''
    uniform = np.random.rand(*args, **kwargs)
    return np.sqrt(2) * sp.erfinv(2 * uniform - 1)
    

if (__name__ == "__main__"):
    np.random.seed(1)

    cor = 0.9

# Problem dimension
    M = 100     # row number of the matrix
    N = 300     # column number

    blkNum = 20 # number of nonzero blocks
    blkLen = 4  # block length

    iterNum = 1

    Phi = randn2(N,M).T

    Phi = Phi / (np.ones((M,1)) * np.sqrt(np.sum(Phi**2, axis=0)))

# Generate nonzero block coefficients
    beta = np.ones((blkNum)) * cor
    blks = np.zeros((blkNum,blkLen))
    blks[:,0] = randn2(blkNum)
    for i in range(1,blkLen):
        blks[:,i] = beta * blks[:,i-1] + np.sqrt(1.0-beta**2) * randn2(blkNum)

    blks = blks / (np.sqrt(np.sum(blks**2, axis=1, keepdims=True)) @ np.ones((1,blkLen)))

    nonzeroCoeff = blks.reshape((blkNum*blkLen,1))


    psbBlockNum = np.floor(N / blkLen).astype('int')

    ind2 = np.array([75,1,55,36,33,54,31, 8,10,61,15,47,27,51,17,32, 2, 4,21,60,43,19,72,29, 9,73, 3,22,20, 6,46,12,49,35,45,28,68,62,67,56,26,13,34,16,18,69,14,70,41,24,64,44,71,58, 7,25,66,52,74, 5,53,11,38,59,42,39,57,63,65,48,37,40,23,50,30])-1
    # ind2 = np.random.permutation(psbBlockNum)

    indice2 = ind2[0:blkNum]
    Wgen = np.zeros((N,1))
    blkLoc = []
    for i in range(blkNum):
        Wgen[ indice2[i] * blkLen : (indice2[i]+1) * blkLen] = nonzeroCoeff[i*blkLen : (i+1)*blkLen]

    Y = Phi @ Wgen

    groupStartLoc = np.arange(0,N,blkLen)

    out = bsbl_b0(Phi, Y, groupStartLoc, 0, verbose=True)

    pl.plot(Wgen)
    pl.plot(out['x'])

    # Observation noise
    # SNR = 15   
    # signal = Phi @ Wgen
    # stdnoise = np.std(signal)*10**(-SNR/20.0)
    # noise = randn2(M,1) * stdnoise

    # Y = signal + noise

    # out = bsbl_b0(Phi, Y, groupStartLoc, 2)