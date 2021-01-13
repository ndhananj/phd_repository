
from autocorrelation import *
from modes import *

def plot_transformed_std(coordinates='transformed_coords.npy'):
    coords = load_matrix(coordinates)
    plot_trans_std(coords)

if __name__ == '__main__':
    coordinates = sys.argv[1] if len(sys.argv)>1 \
        else 'transformed_coords.npy'
    plot_transformed_std(coordinates)
