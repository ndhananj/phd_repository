
import numpy as np
import matplotlib.pyplot as plt

def transformed_coords(coords,S):
    return np.matmul(coords,S.T)

def autocorr1D(x):
    N = x.shape[0]
    A = np.zeros(x.shape[0])
    x -= np.mean(x)
    for k in range(N):
        A[k]=np.sum(np.array([x[i]*x[i+k] for i in range(N-k)]))/N
    return A

def autocorrelate(coords):
    return np.apply_along_axis(autocorr1D, 0, coords)

def autocorrelate_transformed(coords,S):
    return autocorrelate(transformed_coords(coords, S))

def plot_autocorrelate(corr,idx,time_step=0.05):
    N=corr.shape[0]
    T=N*time_step
    subscript = str(idx)+str(idx)
    fig= plt.figure()
    ax = fig.add_subplot(111)
    ax.scatter(np.linspace(0,T,N),corr)
    ax.set_xlabel('time (ps)')
    ax.set_ylabel(r"$K_{"+subscript+"}(\AA^2)$")
    plt.ylim(corr.min()*1.1,corr.max()*1.1)
    plt.savefig("autocorr"+str(idx)+".jpg")
    #plt.show()
