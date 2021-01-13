
from modes import *
from autocorrelation import *

import sys
import glob

def avg_auto_cov(
    globPat='run*/autocorr_', \
    autocorr='avg/autocorr_', \
    start_mode=0, \
    end_mode=9, \
    ):
    for i in range(start_mode,end_mode+1):
        globPat_i = globPat+str(i)+'.npy'
        print(globPat_i)
        filenames = [filename for filename in glob.glob(globPat_i)]
        corr = sum(load_matrix(filename) for filename in filenames)
        corr /= len(filenames)
        print(corr.shape)
        save_matrix(autocorr+str(i)+'.npy', corr)
        plot_autocorrelate(corr,i)

if __name__ == '__main__':
    #inputs
    globPat = sys.argv[1] if len(sys.argv)>1 \
        else '../run*/autocorr_'
    autocorr = sys.arg[2] if len(sys.argv)>2 \
        else 'autocorr_'
    start_mode = sys.argv[3] if len(sys.argv)> 3 \
       else 0
    end_mode = sys.argv[4] if len(sys.argv)>4 \
       else 9
    # calculate
    avg_auto_cov(
        globPat, \
        autocorr, \
        int(start_mode), \
        int(end_mode))
