
from autocorrelation import *
from modes import *

def autocorr(\
    coordinates,\
    mode,\
    time_end=1000, \
    time_start=0,\
    time_step=0.05,\
    ):
    Q = load_matrix(coordinates)
    N=min(int(time_end/time_step),int(Q.shape[0]))
    Qcorr=autocorr1D(Q,N)
    save_matrix('autocorr_'+str(mode)+'.npy',Qcorr[:N])

if __name__ == '__main__':
    coordinates = sys.argv[1]
    mode = int(sys.argv[2])
    time_end = float(sys.argv[3]) if len(sys.argv)>3 else 500
    time_start = float(sys.argv[4]) if len(sys.argv)>4 else 0
    time_step = float(sys.argv[5]) if len(sys.argv)>5 else 0.05
    autocorr(coordinates,mode,time_end,time_start,time_step)
