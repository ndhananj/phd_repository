
import numpy as np
import matplotlib.pyplot as plt

def transformed_coords(coords,S):
    print("coords",coords.shape)
    print("S",S.shape)
    return np.matmul(coords-coords.mean(axis=0),S.T)

def autocorr1D(x,t):
    if t<=0:
        t=N
    print("to_autocorr",x.shape)
    N = x.shape[0]
    A = np.zeros(x.shape[0])
    x -= np.mean(x)
    for k in range(t):
        A[k]=np.mean(np.array([x[i]*x[i+k] for i in range(N-k)]))
    return A/A[0]

def autocorr1D_transformed(coords,S,i):
    return autocorr1D(transformed_coords(coords,S)[:,i])

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
    ax.set_ylabel(r"$K_{"+subscript+"}$")
    plt.ylim(min([0,corr.min()*1.1]),corr.max()*1.1)
    plt.savefig("autocorr"+str(idx)+".jpg")
    #plt.show()

# plot coordinates that are assumed transformed
def plot_transformed(coords,idx,time_step=0.05):
    N=coords.shape[0]
    T=N*time_step
    subscript = str(idx)
    fig= plt.figure()
    ax = fig.add_subplot(111)
    ax.plot(np.linspace(0,T,N),coords,linewidth=0.3,color='k')
    ax.set_xlabel('time (ps)')
    ax.set_ylabel(r"$Q(t)_{"+subscript+"}(\AA)$")
    plt.ylim(min([0,coords.min()*1.1]),coords.max()*1.1)
    plt.savefig("transformed"+subscript+".jpg")

# plot coordinates that are assumed transformed
def plot_transformed_phase(coords,idx,time_step=0.05):
    N=coords.shape[0]
    c=np.reshape(coords,(N))
    vel=(c[1:]-c[:N-1])/time_step
    pos=(c[1:]+c[:N-1])/2
    subscript = str(idx)
    fig= plt.figure()
    ax = fig.add_subplot(111)
    ax.plot(pos,vel,linewidth=0.3,color='k')
    ax.set_xlabel(r"$Q(t)_{"+subscript+"}(\AA)$")
    ax.set_ylabel(r"$\dot{Q}(t)_{"+subscript+"}(\AA/ps)$")
    plt.xlim(min([0,pos.min()*1.1]),pos.max()*1.1)
    plt.ylim(min([0,vel.min()*1.1]),vel.max()*1.1)
    plt.savefig("transformed_phased"+subscript+".jpg")

# plot coordinates that are assumed transformed
def plot_trans_std(coords):
    sigma=coords.std(axis=0)
    N=sigma.shape[0]
    fig= plt.figure()
    ax = fig.add_subplot(111)
    ax.scatter(range(N),sigma)
    ax.set_xlabel(r"$\lambda$");
    ax.set_ylabel(r"$\sigma_{"+"\lambda"+"}(\AA)$")
    plt.ylim(0,coords.max()*1.1)
    plt.savefig("mode_std.jpg")
