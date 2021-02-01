
from autocorrelation import *
from modes import *

def get_autocorrelation(\
    time_start=0,\
    time_step=0.05,\
    coordinates='coordinates.npy',\
    eigenmatrix='eigenmatrix.npy', \
    start_mode=0, end_mode=9\
    ):
    coords = load_matrix(coordinates)
    print('coords',coords.shape)
    S = load_matrix(eigenmatrix)
    idx_start=int(time_start/time_step)
    means=coords[idx_start:,:].mean(axis=0)
    print('means',mean.shape)
    Q=transformed_coords(coords-means,S)
    save_matrix('transformed_coords.npy', Q)
    N=min(40000,int(Q.shape[0]/5))
    for i in range(int(start_mode),int(end_mode)+1):
        Qcorr=autocorr1D(Q[:,i],N)
        save_matrix('autocorr_'+str(i)+'.npy',Qcorr[:N])
        plot_autocorrelate(Qcorr[:N],i)

if __name__ == '__main__':
    time_start = float(sys.argv[1]) if len(sys.argv)>1 else 0
    time_step = float(sys.argv[2]) if len(sys.argv)>2 else 0.05
    coordinates = sys.argv[3] if len(sys.argv)>3 else 'coordinates.npy'
    eigenmatrix = sys.argv[4] if len(sys.argv)>4 else 'eigenmatrix.npy'
    start_mode = sys.argv[5] if len(sys.argv)>5 else 0
    end_mode = sys.argv[6] if len(sys.argv)>6 else 9
    get_autocorrelation(\
        time_start,time_step,coordinates,eigenmatrix,start_mode,end_mode)
