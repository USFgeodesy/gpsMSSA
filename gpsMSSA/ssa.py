'''
Singular Spectrum Analysis for a single component of GPS station

'''
import numpy as np
import matplotlib.pyplot as plt

from scipy.linalg import toeplitz

from oceans.ff_tools import lagcorr

def shift(arr, n, order='forward'):
    '''
    function for creating lagged array from original Time series
    '''
    if order == 'forward':
        shifted = arr[n:] + [0] * n
    elif order == 'reversed':
        shifted = [0] * n + arr[:-n]
    else:
        print("Order %s not recognized.  Try forward or reversed" % order)

    return shifted

def embedTS(ts,M):
    '''
    create embedded time series matrix Y
    inputs:
        ts = timeseries
        M = window length
    otuputs
        Y = embedded matrix
    '''
    Y = list(ts)
    for i in np.arange(1,M):
        Y = np.c_[Y,shift(list(ts),i)]
    return Y


def eigvalues(ts,M,plot=False):
    '''
    get the eigenvalues from a single time series
    ts = GPS time series
    M = window length
    plot = Boolean (default = False) plot the eigvalues
    returns:
        eigenvalues, eigenvectors
    '''
    covX = lagcorr(ts, ts, M=M)
    C = toeplitz(covX)
    lamb, rho = np.linalg.eig(C)
    if plot==True:
        fig, ax = plt.subplots(figsize=(8, 6))
        ax.plot(range(M), rho, marker='.')
        ax.legend(['1', '2', '3', '4'], numpoints=1)
        ax.set_title(r'Eigenvectors "$\rho$"')
        ax.plot(range(M), [0] * len(rho), 'k--')
        plt.show()
    return lamb,rho

def pComp(Y,rho):
    '''
    calculate prinicpal component array
    inputs:
        Y : embedded lagged time series
        rho: Eigenvectors
    return:
        pc = array of principal components
    '''
    pc = np.dot(Y, rho)
    return pc

def plotPC(pc,ts,pcArray):
    '''
    plot the principal components in pcArray
    inputs:
        pc = array of princapal components
        ts = timeseries
        pcArray = indices of principal components to plot
    '''
    # Plotting the PCs.
    fig, ax = plt.subplots(nrows=len(pcArray), ncols=1, figsize=(8, 6))
    ax[0].set_title(r'Principal components PC vs T')
    for i in range(len(pcArray)):
        ax[i].plot(ts.index, pc[:,i], '.-')
    return

def reconstruct_ts(PC,rho,indexes,M, plot = False):
    '''
    reconstruct the original time series
    inputs:
    PC =  array of principal components
    rho = eigenvector array
    indexes = indexes of principal comonents to recreate
    M = window length
    optional :
        plot : plot the reconstructed time series (default = False)
    outputs:
        array of reconstructed time series from indexes references
    '''
    Zarray = []
    for i in indexes:
        Z = list(PC[:,i])
        for ind in np.arange(1,M):
            Z = np.c_[Z,shift(list(PC[:,i]),ind,'reversed')]
        Zarray.append(Z)
    RCarray = []
    for i,value in enumerate(indexes):
        RCarray.append(np.dot(Zarray[i], rho[:, value]) / M)
    if plot==True:
        fig, ax = plt.subplots(nrows=len(indexes), ncols=1, figsize=(8, 6))
        ax[0].set_title('Reconstructed components RC vs T')
        for i in range(len(indexes)):
            ax[i].plot(range(PC[:,i]), RCarray[i], 'r.-',label = str(indexes[i]))
        plt.show()
    return RCarray

def compareRC_orig(ts,RCList):
    '''
    Compare linear combination of components with original time series
    ts, list of components to compare
    '''
    RCtot = np.zeros(len(RCList[0]))
    for i in range(len(RCList)):
        RCtot = RCtot + RCList[i]
    fig, (ax0) = plt.subplots(nrows=1, ncols=1, figsize=(8, 6))
    ax0.set_title(r'Reconstructed components RC vs T')
    ax0.plot(ts.index, ts.get_values(), 'b.-', label='Original')
    ax0.plot(ts.index, RCtot, 'r.-', label='RC')
    ax0.legend(numpoints=1)
    plt.show()

def ssa_app(ts,M,prinComp):
    # generate Y matrix
    Y = embedTS(ts,M)
    # get eigenvectors and eigenvalues
    lamb,rho = eigvalues(ts,M)
    #get pC
    pc = pComp(Y,rho)
    #reconstruct Time Series
    rc = reconstruct_ts(pc,rho,prinComp,M)
    #Compare
    compareRC_orig(ts,rc)
    return
