
from autocorrelation import *
from modes import *

def get_transformed(\
    time_start=0,\
    time_step=0.05,\
    coordinates='coordinates.npy',\
    eigenmatrix='eigenmatrix.npy', \
    ):
    coords = load_matrix(coordinates)
    S = load_matrix(eigenmatrix)
    idx_start=int(time_start/time_step)
    means=coords[idx_start:,:].mean(axis=0).reshape((1,coords.shape[1]))
    print('means',means.shape)
    Q=transformed_coords(coords-means,S)
    save_matrix('transformed_coords.npy', Q)

if __name__ == '__main__':
    time_start = float(sys.argv[1]) if len(sys.argv)>1 else 0
    time_step = float(sys.argv[2]) if len(sys.argv)>2 else 0.05
    coordinates = sys.argv[3] if len(sys.argv)>3 else 'coordinates.npy'
    eigenmatrix = sys.argv[4] if len(sys.argv)>4 else 'eigenmatrix.npy'
    get_transformed(time_start,time_step,coordinates,eigenmatrix)
