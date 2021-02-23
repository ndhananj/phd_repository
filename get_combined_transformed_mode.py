
from autocorrelation import *
from modes import *
import os.path
from os import path
import tables

def get_transformed_mode(coordinates='transformed_coords.npy',mode=0):
    i=0
    names = []
    while path.exists(chunk_name(coordinates,i)):
        names.append(chunk_name(coordinates,i))
        i=i+1
    coords = (load_matrix(name) for name in names)
    used = (c[:,mode] for c in coords)
    save_matrix("Q"+str(mode)+"(t).npy",np.concatenate(list(used)))

if __name__ == '__main__':
    coordinates = sys.argv[1] if len(sys.argv)>1 else 'transformed_coords.npy'
    mode = int(sys.argv[2]) if len(sys.argv)>2 else 0
    get_transformed_mode(coordinates,mode)
