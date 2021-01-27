
from autocorrelation import *
from modes import *
import os.path
from os import path

def plot_autocorr(\
    coordinates='transformed_coords.npy',\
    start_mode=0, end_mode=9,time_end=0,time_step=0.05
    ):
    i=0
    names = []
    while path.exists(chunk_name(coordinates,i)):
        names.append(chunk_name(coordinates,i))
        i=i+1
    for j in range(start_mode,end_mode+1):
        coords = (load_matrix(name) for name in names)
        used = (c[:,j] for c in coords)
        corr=chunked_autocorr1D(used,int(time_end/time_step))
        save_matrix('autocorr_'+str(j)+'.npy',corr)
        plot_autocorrelate(corr,j)

if __name__ == '__main__':
    coordinates = sys.argv[1] if len(sys.argv)>1 \
        else 'transformed_coords.npy'
    start_mode = int(sys.argv[2]) if len(sys.argv)>2 \
       else 0
    end_mode = int(sys.argv[3]) if len(sys.argv)>3 \
       else 9
    time_end = float(sys.argv[4]) if len(sys.argv)>4 \
        else 200
    time_step = float(sys.argv[5]) if len(sys.argv)>5 \
        else 0.05
    plot_autocorr(coordinates,start_mode,end_mode,time_end)
