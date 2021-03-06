
from autocorrelation import *
from modes import *

def plot_transformed_coords(\
    coordinates='transformed_coords.npy',\
    start_mode=0, end_mode=9,time_start=0,time_step=0.05):
    coords = load_matrix(coordinates)
    for i in range(start_mode,end_mode+1):
        plot_transformed(coords[:,i],i,time_step=time_step,time_start=time_start)

if __name__ == '__main__':
    coordinates = sys.argv[1] if len(sys.argv)>1 \
        else 'transformed_coords.npy'
    start_mode = int(sys.argv[2]) if len(sys.argv)>2 \
       else 0
    end_mode = int(sys.argv[3]) if len(sys.argv)>3 \
       else 9
    time_start = float(sys.argv[4]) if len(sys.argv)>4 \
        else 0
    time_step = float(sys.argv[5]) if len(sys.argv)>5 \
        else 0.05
    plot_transformed_coords(coordinates,start_mode,end_mode,time_start,time_step)
