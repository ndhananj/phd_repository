
from autocorrelation import *
from modes import *

def plot_transformed_coords(\
    coordinates='transformed_coords.npy',\
    start_mode=0, end_mode=9,time_start=1500):
    coords = load_matrix(coordinates)
    for i in range(start_mode,end_mode+1):
        plot_transformed(coords[:,i],i)

if __name__ == '__main__':
    coordinates = sys.argv[1] if len(sys.argv)>1 \
        else 'transformed_coords.npy'
    start_mode = sys.argv[2] if len(sys.argv)>2 \
       else 0
    end_mode = sys.argv[3] if len(sys.argv)>3 \
       else 9
    time_start = float(sys.argv[4]) if len(sys.argv)>3 \
        else 1500
    plot_transformed_coords(coordinates,start_mode,end_mode,time_start)
