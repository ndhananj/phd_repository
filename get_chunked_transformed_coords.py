
from autocorrelation import *
from modes import *
import os.path
from os import path

def get_transformed(\
    time_start=0,\
    time_step=0.05,\
    coordinates='coordinates.npy',\
    eigenmatrix='eigenmatrix.npy', \
    means='means.npy'
    ):
    S = load_matrix(eigenmatrix)
    mean = load_matrix(means)
    mean = mean.reshape((1,mean.shape[0]))
    print('mean',mean.shape)
    i=0
    while path.exists(chunk_name(coordinates,i)):
        coords = load_matrix(chunk_name(coordinates,i))
        if 1==len(coords.shape):
            coords=coords.reshape(1,coords.shape[0])
        idx_start=int(time_start/time_step)
        Q=transformed_coords(coords-mean,S)
        save_matrix(chunk_name('transformed_coords.npy',i), Q)
        i=i+1

if __name__ == '__main__':
    time_start = float(sys.argv[1]) if len(sys.argv)>1 else 0
    time_step = float(sys.argv[2]) if len(sys.argv)>2 else 0.05
    coordinates = sys.argv[3] if len(sys.argv)>3 else 'coordinates.npy'
    eigenmatrix = sys.argv[4] if len(sys.argv)>4 else 'eigenmatrix.npy'
    means = sys.argv[5] if len(sys.argv)>5 else 'means.npy'
    get_transformed(time_start,time_step,coordinates,eigenmatrix,means)
