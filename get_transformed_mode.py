
from autocorrelation import *
from modes import *
import os.path
from os import path
import tables

def get_transformed_mode(coordinates='transformed_coords.npy',mode=0):
    coords=load_matrix(coordinates)
    save_matrix("Q"+str(mode)+"(t).npy",coords[:,mode])

if __name__ == '__main__':
    coordinates = sys.argv[1] if len(sys.argv)>1 else 'transformed_coords.npy'
    mode = int(sys.argv[2]) if len(sys.argv)>2 else 0
    get_transformed_mode(coordinates,mode)
