
from autocorrelation import *
from modes import *

def get_autocorrelation(\
    time_end=1000, \
    time_start=0,\
    time_step=0.05,\
    coordinates='transformed_coords.npy',\
    start_mode=0, end_mode=9\
    ):
    Q = load_matrix(coordinates)
    N=min(int(time_end/time_step),int(Q.shape[0]))
    for i in range(int(start_mode),int(end_mode)+1):
        Qcorr=autocorr1D(Q[:,i],N)
        save_matrix('autocorr_'+str(i)+'.npy',Qcorr[:N])
        #plot_autocorrelate(Qcorr[:N],i)

if __name__ == '__main__':
    time_end = float(sys.argv[1]) if len(sys.argv)>1 else 500
    time_start = float(sys.argv[2]) if len(sys.argv)>2 else 0
    time_step = float(sys.argv[3]) if len(sys.argv)>3 else 0.05
    coordinates = sys.argv[4] if len(sys.argv)>4 else 'transformed_coords.npy'
    start_mode = sys.argv[5] if len(sys.argv)>5 else 0
    end_mode = sys.argv[6] if len(sys.argv)>6 else 9
    get_autocorrelation(\
        time_end,time_start,time_step,coordinates,start_mode,end_mode)
