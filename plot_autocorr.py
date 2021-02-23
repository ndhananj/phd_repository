
from autocorrelation import *
from modes import *

def plot_autocorr(\
    corr='autocorr_',\
    start_mode=0, end_mode=9,time_start=0):
    for i in range(start_mode,end_mode+1):
        Qcorr = load_matrix(corr+str(i)+'.npy')
        plot_autocorrelate(Qcorr,i)

if __name__ == '__main__':
    corr = sys.argv[1] if len(sys.argv)>1 \
        else 'autocorr_'
    start_mode = int(sys.argv[2]) if len(sys.argv)>2 \
       else 0
    end_mode = int(sys.argv[3]) if len(sys.argv)>3 \
       else 9
    time_start = float(sys.argv[4]) if len(sys.argv)>4 \
        else 0
    plot_autocorr(corr,start_mode,end_mode,time_start)
