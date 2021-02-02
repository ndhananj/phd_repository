
import numpy as np
import matplotlib.pyplot as plt
import scipy
import scipy.optimize
import vaex

# assumes that the coordinated have alrady been set to a mean
def transformed_coords(coords,S):
    print("coords",coords.shape)
    print("S",S.shape)
    return np.matmul(coords,S)

def corr_comps(x,k,N):
    return (x[i]*x[i+k] for i in range(N-k))

def autocorr1D(x,t):
    print("to_autocorr",x.shape)
    N = x.shape[0]
    if not(t) or t<=0:
        t=N
    A = np.zeros(t)
    x -= np.mean(x)
    for k in range(t):
        A[k]=np.mean(np.array(list(corr_comps(x,k,N))))
    return A/A[0]

def chunked_autocorr1D(chunks,t):
    return autocorr1D(np.concatenate(list(chunks)),t)

def autocorr1D_transformed(coords,S,i):
    return autocorr1D(transformed_coords(coords,S)[:,i])

def autocorrelate(coords):
    return np.apply_along_axis(autocorr1D, 0, coords)

def autocorrelate_transformed(coords,S):
    return autocorrelate(transformed_coords(coords, S))

def sum_simple_exp_decay(t,a,tau1,tau2):
    return a*np.exp(-t/tau1)+(1-a)*np.exp(-t/tau2)

def format_label(decimal):
    return "{:.2f}".format(decimal)

def sum_simple_exp_label(a,tau1,tau2):
    return "$"+format_label(a)+"e^{-t/"+format_label(tau1)+"}+"\
     +format_label(1-a)+"e^{-t/"+format_label(tau2)+"}$"#+format_label(c)+"$"

def fit_sum_exp(corr,time,f=sum_simple_exp_decay):
    popt, pcov = scipy.optimize.curve_fit(f,time,corr)
    return popt, pcov

def plot_autocorrelate(corr,idx,time_step=0.05):
    N=corr.shape[0]
    T=N*time_step
    subscript = str(idx)+str(idx)
    time = np.linspace(0,T,N)
    popt, pcov = fit_sum_exp(corr,time,f=sum_simple_exp_decay)
    fig= plt.figure()
    ax = fig.add_subplot(111)
    ax.plot(time,corr,linewidth=0.3,color='k',label='values')
    ax.plot(
       time,sum_simple_exp_decay(time,*popt),
       linewidth=0.3,color='k',linestyle='--',
       label=sum_simple_exp_label(*popt)
    )
    ax.set_xlabel('time (ps)')
    ax.set_ylabel(r"$K_{"+subscript+"}$")
    plt.ylim(min([0,corr.min()*1.1]),corr.max()*1.1)
    plt.legend()
    plt.savefig("autocorr"+str(idx)+".png",dpi=600)
    #plt.show()


# plot coordinates that are assumed transformed
def plot_transformed(
    coords,idx,time_step=0.05,time_start=0,basename="transformed"
    ):
    print(coords)
    N=coords.shape[0]
    T=N*time_step
    w=np.convolve(coords, np.ones(2000)/2000, mode='same')
    times=np.linspace(0,T,N)
    subscript = str(idx)
    fig= plt.figure()
    ax = fig.add_subplot(111)
    ax.plot(times,coords,linewidth=0.3,color='b',linestyle='--',alpha=0.5)
    ax.plot(times,w,linewidth=1,color='r')
    ax.set_xlabel('time (ps)')
    ax.set_ylabel(r"$Q(t)_{"+subscript+"}(\AA)$")
    plt.ylim(min([0,coords.min()*1.1]),coords.max()*1.1)
    #plt.show()
    plt.savefig(basename+subscript+".png",dpi=600)

# plot coordinates that are assumed transformed
def plot_chunked_transformed(
    chunked_coords,idx,time_step=0.05,time_start=0,basename="transformed"
    ):
    subscript = str(idx)
    fig= plt.figure()
    ax = fig.add_subplot(111)
    ymin=0
    ymax=0
    for coords in chunked_coords:
        if(len(coords)>1):
            print(coords)
            N=coords.shape[0]
            T=N*time_step
            w=np.convolve(coords, np.ones(2000)/2000, mode='same')
            times=np.linspace(0,T,N)+time_start
            ax.plot(times,coords,linewidth=0.3,color='b',linestyle='--',alpha=0.5)
            ax.plot(times,w,linewidth=1,color='r')
            time_start+=T
        else:
            print("Short coords")
        ymin = min([ymin,coords.min()*1.1])
        ymax = max([ymax,coords.max()*1.1])
    ax.set_xlabel('time (ps)')
    ax.set_ylabel(r"$Q(t)_{"+subscript+"}(\AA)$")
    plt.ylim(ymin,ymax)
    #plt.show()
    plt.savefig(basename+subscript+".png",dpi=600)

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
    #plt.show()
    plt.savefig("transformed_phased"+subscript+".jpg")

# plot coordinates that are assumed transformed
def plot_chunked_transformed_phase(
    chunked_coords,idx,time_step=0.05,time_start=0,basename="transformed_phase"
    ):
    subscript = str(idx)
    fig= plt.figure()
    ax = fig.add_subplot(111)
    ymin, ymax, xmin, xmax = 0, 0, 0, 0
    for coords in chunked_coords:
        if(len(coords)>1):
            print(coords)
            N=coords.shape[0]
            c=np.reshape(coords,(N))
            vel=(c[1:]-c[:N-1])/time_step
            pos=(c[1:]+c[:N-1])/2
            ax.plot(pos,vel,linewidth=0.3,color='k')
            xmin = min([xmin,pos.min()*1.1])
            xmax = max([xmax,pos.max()*1.1])
            ymin = min([ymin,vel.min()*1.1])
            ymax = max([ymax,vel.max()*1.1])
    ax.set_xlabel(r"$Q(t)_{"+subscript+"}(\AA)$")
    ax.set_ylabel(r"$\dot{Q}(t)_{"+subscript+"}(\AA/ps)$")
    plt.xlim(xmin,xmax)
    plt.ylim(ymin,ymax)
    #plt.show()
    plt.savefig(basename+subscript+".png",dpi=600)

# plot coordinates that are assumed transformed
def plot_chunked_transformed_phase_hist(
    chunked_coords,idx,basename="transformed_phase_hist", time_step=0.05
    ):
    subscript = str(idx)
    fig= plt.figure()
    ax = fig.add_subplot(111)
    ymin, ymax, xmin, xmax = 0, 0, 0, 0
    vel_list = []
    pos_list = []
    for coords in chunked_coords:
        if(len(coords)>1):
            print(coords)
            N=coords.shape[0]
            c=np.reshape(coords,(N))
            vel = (c[1:]-c[:N-1])/time_step
            pos = (c[1:]+c[:N-1])/2
            vel_list.extend(vel)
            pos_list.extend(pos)
            #ax.plot(pos,vel,linewidth=0.3,color='k')
            xmin = min([xmin,pos.min()*1.1])
            xmax = max([xmax,pos.max()*1.1])
            ymin = min([ymin,vel.min()*1.1])
            ymax = max([ymax,vel.max()*1.1])
    h = ax.hist2d(pos_list,vel_list,bins=40)
    fig.colorbar(h[3], ax=ax)
    ax.set_xlabel(r"$Q(t)_{"+subscript+"}(\AA)$")
    ax.set_ylabel(r"$\dot{Q}(t)_{"+subscript+"}(\AA/ps)$")
    plt.xlim(xmin,xmax)
    plt.ylim(ymin,ymax)
    #plt.show()
    plt.savefig(basename+subscript+".png",dpi=600)

# plot coordinates that are assumed transformed
def plot_chunked_transformed_hist(
    chunked_coords,idx,basename="transformed_hist", time_step=0.05
    ):
    subscript = str(idx)
    fig= plt.figure()
    ax = fig.add_subplot(111)
    ymin, ymax, xmin, xmax = 0, 0, 0, 0
    pos_list = []
    for coords in chunked_coords:
        if(len(coords)>1):
            print(coords)
            N=coords.shape[0]
            c=np.reshape(coords,(N))
            pos = c
            pos_list.extend(pos)
            #ax.plot(pos,vel,linewidth=0.3,color='k')
            xmin = min([xmin,pos.min()*1.1])
            xmax = max([xmax,pos.max()*1.1])
    h = ax.hist(pos_list,bins=40)
    #fig.colorbar(h)
    ax.set_xlabel(r"$Q(t)_{"+subscript+"}(\AA)$")
    ax.set_ylabel("count")
    plt.xlim(xmin,xmax)
    #plt.show()
    plt.savefig(basename+subscript+".png",dpi=600)

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
