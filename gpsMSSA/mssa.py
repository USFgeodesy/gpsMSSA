'''
Multi-variate single spectrum analysis for GPS network
'''

import numpy as np
import matplotlib.pylab as plt
from scipy.linalg import toeplitz
from oceans.ff_tools import lagcorr
import gpsUtils as gps
import ssa
import glob
from scipy.signal import lfilter, filtfilt

def create_network_matrix(network_folder, start, end):
    '''
    parse network folder to make network matrix
    GPS1Lat GPS1Lon GPS1Vert GPS2Lat ...
    all time series will be same length and be filled, and zscored
    '''

    #loop through directory
    fs = glob.glob(network_folder + '/*')
    for i, f in enumerate(fs):
        print f
        ts = gps.readData(f)
        sr = ts.ts2pandas(start, end)
        if i == 0:
            df = sr.to_frame(name =f.split('/')[1])
        else:
            df[f.split('/')[1]] = sr
    return df

def create_Y(df, M):
    '''
        create embedded time series matrix Y
        inputs:
            X = data frame with data
            M = window length
        otuputs
            Y = embedded matrix
    '''
    for i, c in enumerate(df.columns):
        R = list(df[c])
        for u in np.arange(1, M):
            R = np.c_[R, ssa.shift(list(df[c]), u)]
        if i == 0:
            Y = R
        else:
            Y = np.c_[Y, R]
    C = np.dot(Y.T, Y) / len(df)
    return Y, C

def toep_cov(df, M):
    '''
    compute the Block Toeplitz covariance method instead of using lagged copies

    '''
    # ctrate the covariance matrix
    N = len(df)
    L = len(df.columns)
    X = df.as_matrix()
    c = np.zeros((2*M+1, L**2))
    #nfft=2.^(fix(log(N-1)/log(2))+2); % twice the next power of two
    nfft = int(2**(np.fix(np.log(N-1)/np.log(2))+2))
    xtmp = np.zeros((nfft,L))
    for i in range(L):
        xtmp[:,i]= np.fft.fft(X[:,i],nfft,axis = 0 )
    ind = np.arange(nfft/2-M,nfft/2+M+1)
    for i in range(L):
        for j in range(i,L):
            y = np.fft.ifft(xtmp[:,j]*np.conj(xtmp[:,i]),axis = 0 )
            y = np.real(np.fft.fftshift(y))
            col = (i)*L+j
            c[:,col]=y[ind]
    ind = np.arange(2*M+1,0,-1)
    for i in range(1,L):
        for j in range(i-1):
            col = (i-1)*L+j
            colnew = (j)*L+i
            c[:,col]=c[ind-1,colnew]
    nrm = np.concatenate((np.arange(N-M,N), np.arange(N,N-M-1,-1)))
    c = c/np.tile(nrm,(L**2,1)).T
    # build the T covariance matroc bloc by Block
    C = np.zeros((L*M,M*L))
    #diagonal blocks
    for i in range(L):
        q = c[M-1:2*M-1,(i)*(L+1)]
        Tij = toeplitz(q)
        C[i*M:(i+1)*M,(i)*M:(i+1)*M]= Tij
    #build all the other blocks
    for i in range(L):
        for j in range(i+1,L):
            q=c[:,(i)*L+j];
            Tij=toeplitz(q[M:0:-1],q[M-1:2*M-1]);
            print np.shape(Tij)
            C[(i)*M:(i+1)*M, (j)*M:(j+1)*M]=Tij  # Above diagonal
            C[(j)*M:(j+1)*M, (i)*M:(i+1)*M]=Tij.T # Below diagonal
    C = np.dot(C,C.T)/2.0
    return C



def get_eigen(C, plot = False):
    '''
    return eigenvalues and eigenvectors of C
    '''
    lamb, rho = np.linalg.eig(C)
    idx = lamb.argsort()[::-1]
    rho = rho[:, idx]
    lamb = lamb[idx]
    #rho[:, 0] = -rho[:, 0]
    if plot==True:
        fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(8, 6))
        ax.set_title(r'Eigenvectors "$\rho$"')
        ax.plot(rho[:, 0], 'b.-', label='1')
        ax.plot(rho[:, 1], 'g.-', label='2')
        ax.legend(numpoints=1)
        plt.show()
    return rho, lamb

def get_PC(rho,M,df,plot=False):
    '''
    get principal components of matrix
    '''
    E = rho
    X = df.as_matrix()
    L = len(df.columns)
    A = np.zeros((len(df),np.shape(rho)[1]))
    a = np.zeros((len(df),len(df.columns)))
    Ej = np.zeros((M,len(df.columns)))
    for K in range(np.shape(E)[1]):
        Ej = Ej.flatten()
        Ej = E[:,K]
        Ej = Ej.reshape((M,len(df.columns)))
        Ej = np.flipud(Ej)
        for j in range(L):
            a[:,j]= lfilter(Ej[:,j],1,X[:,j],axis = 0 )
            #a[:,j]= filtfilt(Ej[:,j],1,X[:,j],axis = 0,method = 'pad')
        if L>1:
            A[:,K] = (np.sum(a.T,axis = 0))
        else:
            A[:,K] = a
    A = A[M-1:len(df),:]

    if plot==True:
        fig, (ax0, ax1, ax2, ax3) = plt.subplots(nrows=4, ncols=1, figsize=(8, 6))
        ax0.set_title(r'Principal components PC vs t')
        ax0.plot(PC[:, 0], 'b.-')
        ax1.plot(PC[:, 1], 'b.-')
        ax2.plot(PC[:, 2], 'b.-')
        ax3.plot(PC[:, 3], 'b.-')
        plt.show()

    return A

def plotPC(pc,ts,pcArray):
    '''
    plot the principal components in pcArray
    inputs:
        pc = array of princapal components
        ts = timeseries
        pcArray = indices of principal components to plot
    '''
    # Plotting the PCs.
    fig, ax = plt.subplots(nrows=len(pcArray), ncols=1, figsize=(15, 15))
    ax[0].set_title(r'Principal components PC vs T')
    for i,ix in enumerate(pcArray):
        ax[i].plot(ts.index, pc[:,ix], '.-')
        ax[i].set_ylabel('PC ' + str(pcArray[i]))
        #ax[i].axvspan('03-07-11','06-14-11',alpha = 0.5)
    return


def recreat_ts(df,PC,rho,M,station):
    '''
    recreate the time series by reprojecting
    PC back on the eigenvectors
    '''
    # Reconstruction of the time series.
    N = len(df)
    X0 = df[station]
    RC = np.zeros((N, N-M+1))
    index = df.columns.get_loc(station)
    for k in np.arange(M):
        print k
        for i in np.arange(1,M):
                xi = 0
                for j in np.arange(i):
                    xi=xi + PC[k,i-j+1]*rho[j,k]
                RC[i,k]=xi/i
        for i in np.arange(M,N-M+1):
                xi = 0
                for j in np.arange(M):
                    xi = xi + PC[k,i-j+1]*rho[j,k]
                RC[i,k] = xi/M
        for i in np.arange(N-M+2,N):
            xi = 0
            for j in np.arange(i-N+M,M):
                xi = xi + PC[k,i-j+1]*PC[j,k]
            RC[i,k] = xi/(N-i+1)
    RC0 = RC
    '''
    for m in np.arange(M):
        Z = np.zeros((N, M))  # Time-delayed embedding of PC[:, m].
        for m2 in np.arange(M):
            Z[m2 - N:, m2] = PC[:N - m2, m]

        # Determine RC as a scalar product.
        RC0[:, m] = np.dot(Z, rho[index*M:(index+1)*M, m] / M
    '''
    fig, axs = plt.subplots(nrows=4, ncols=1, figsize=(8, 6))
    fig.suptitle(r'Reconstruction components RC vs t')
    axs[0].plot(df.index,RC0[:, 0], 'r.-')
    axs[1].plot(df.index,RC0[:, 1], 'r.-')
    axs[2].plot(df.index,RC0[:, 2], 'r.-')
    axs[3].plot(df.index,RC0[:, 3], 'r.-')
    #plt.show()
    plt.savefig(station + '.rcvt.png')
    plt.close()

    """The first both RC1 and RC2 describe the oscillations, where the RC3 and RC4
    describe a trend (which may be introduced due to the random number generator
    and the very short time series).  When we summarize all eight RCs we
    reconstruct the whole time series."""

    fig, axs = plt.subplots(nrows=1, ncols=2, figsize=(8, 6))
    fig.suptitle(r'Original time series and reconstructions vs t')
    axs[0].plot(df.index,X0, '*', label='X0 Original')
    axs[0].plot(df.index,RC0.sum(axis=1), 'r.-', label='RC0-7')
    axs[0].legend(numpoints=1)


    axs[1].plot(df.index,X0, '*', label='X0 Original')
    axs[1].plot(df.index,RC0[:, 0] +RC0[:, 1] + RC0[:, 3], 'r.-', label='RC0-1')
    axs[1].legend(numpoints=1)
    plt.savefig(station + '.orvr.png')
    plt.close()
    #plt.show()
    return RC0

def compareRC_orig(df,station,rc,RCList):
    '''
    Compare linear combination of components with original time series
    ts, list of components to compare
    '''
    RCtot = np.zeros(len(rc))
    for i in RCList:
        RCtot = RCtot + rc[:,i]
    fig, (ax0) = plt.subplots(nrows=1, ncols=1, figsize=(10, 10))
    ax0.set_title(r'Reconstructed components RC vs T')
    ax0.plot(df.index, df[station].get_values(), '*', label='Original')
    ax0.plot(df.index, RCtot, 'r.-', label='RC')
    ax0.legend(numpoints=1)
    #ax0.axvspan('03-07-11','06-14-11',alpha = 0.5)
    plt.savefig(station  + 'rc.png')
    plt.close()
    #plt.show()
    return RCtot
def rc(E,A,L):
    ml,k=np.shape(E)
    ra, ka=np.shape(A)

    M=ml/L;  # %     These lines assume that E and A
    N=ra+M-1 # %     have the "right" row dimensions.

    R=np.zeros((N,L*k));
    Z=np.zeros((M-1,k));
    A=np.concatenate((A.T,Z.T),axis = 1) # % zero pad A to make it N by k
    A=A.T #;
    Ej=np.zeros((M,L))

    #% Calculate RCs
    for j in range(k):
        Ej=E[:,j] # % Convert the j-th EOF into a matrix of filters
        Ej = Ej.reshape((M,L))
        for i in range(L):    # % Compute the RCs for the j-th EOF/PC pair
            R[:,j*L+i]=lfilter(Ej[:,i],M,A[:,j]);
            #R[:,j*L+i]=filtfilt(Ej[:,i],M,A[:,j],method = 'pad');
    #% Adjust first M-1 rows and last M-1 rows
    for i in range(M-1):
      R[i,:]=R[i,:]*(M/(i+1.0))
      R[N-1-(i),:]=R[N-1-(i),:]*M/(i+1)
    Rsum = np.zeros((N,L))
    for j in range(k):
                Rsum = Rsum+R[:,(j)*L:(j+1)*L]
    return R,Rsum
def mc_test(minEig):
    return
