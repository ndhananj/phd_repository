
from autocorrelation import *
from modes import *
import os.path
from os import path

def plot_transformed_coords(\
    coordinates='transformed_coords.npy',\
    start_mode=0, end_mode=9,time_start=0):
    i=0
    while path.exists(chunk_name(coordinates,i)):
        coords = load_matrix(chunk_name(coordinates,i))
        for j in range(start_mode,end_mode+1):
            plot_transformed(coords[:,j],j,basename=chunk_name("transformed",i))
        i=i+1

if __name__ == '__main__':
    coordinates = sys.argv[1] if len(sys.argv)>1 \
        else 'transformed_coords.npy'
    start_mode = int(sys.argv[2]) if len(sys.argv)>2 \
       else 0
    end_mode = int(sys.argv[3]) if len(sys.argv)>3 \
       else 9
    time_start = float(sys.argv[4]) if len(sys.argv)>4 \
        else 0
    plot_transformed_coords(coordinates,start_mode,end_mode,time_start)
