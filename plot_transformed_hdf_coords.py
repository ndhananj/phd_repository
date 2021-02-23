
from autocorrelation import *
from modes import *

def plot_transformed_coords(\
    coordinates='transformed_coords',\
    start_mode=0, end_mode=9,time_start=0):
    coords = load_hdf(coordinates)
    needed_coords=(c[start_mode:end_mode+1] for c in coords)
    for i in range(end_mode+1-start_mode):
        my_coords = np.array([c[i] for c in needed_coords])
        plot_transformed(my_coords,i+start_mode)

if __name__ == '__main__':
    coordinates = sys.argv[1] if len(sys.argv)>1 \
        else 'transformed_coords'
    start_mode = int(sys.argv[2]) if len(sys.argv)>2 \
       else 0
    end_mode = int(sys.argv[3]) if len(sys.argv)>3 \
       else 9
    time_start = float(sys.argv[4]) if len(sys.argv)>4 \
        else 0
    plot_transformed_coords(coordinates,start_mode,end_mode,time_start)
