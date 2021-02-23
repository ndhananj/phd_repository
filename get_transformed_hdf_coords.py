
from autocorrelation import *
from modes import *

def get_transformed(\
    coordinates='coordinates',\
    means='means.npy',\
    eigenmatrix='eigenmatrix.npy', \
    ):
    coords = load_hdf(coordinates)
    S = load_matrix(eigenmatrix)
    my_means=load_matrix(means)
    my_means=my_means.reshape((1,my_means.shape[0]))
    print('means',my_means.shape)
    residues = (c-my_means for c in coords)
    Q=(transformed_coords(r,S) for r in residues)
    save_hdf('transformed_coords', Q)

if __name__ == '__main__':
    coordinates = sys.argv[1] if len(sys.argv)>1 else 'coordinates'
    means = sys.argv[2] if len(sys.argv)>2 else 'means.npy'
    eigenmatrix = sys.argv[3] if len(sys.argv)>3 else 'eigenmatrix.npy'
    get_transformed(coordinates,means,eigenmatrix)
