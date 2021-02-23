
from autocorrelation import *
from modes import *

def plot_Q(\
    coordinates,\
    mode=0,time_start=0,time_step=0.05):
    coords = load_matrix(coordinates)
    plot_transformed(coords,mode,time_step=time_step,time_start=time_start)

if __name__ == '__main__':
    coordinates = sys.argv[1]
    mode = int(sys.argv[2]) if len(sys.argv)>2 \
       else 0
    time_start = float(sys.argv[3]) if len(sys.argv)>3 \
        else 0
    time_step = float(sys.argv[4]) if len(sys.argv)>4 \
        else 0.05
    plot_Q(coordinates,mode,time_start,time_step)
