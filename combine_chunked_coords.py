from autocorrelation import *
from modes import *
import os.path
from os import path
import tables as tb

def combine_coords(coordinates='coordinates.npy',combined='coordinates.dbm'):
    i=0
    f= tb.open_file(combined, mode='w')
    coord = tb.Float32Col()
    while path.exists(chunk_name(coordinates,i)):
        print('processing ',chunk_name(coordinates,i))
        arr=load_matrix(chunk_name(coordinates,i))
        if 0==i :
            ROW_SIZE = arr.shape[1]
            farr = f.create_earray(f.root, 'coordinates', coord, (0, ROW_SIZE))
        for j in range(arr.shape[0]):
            farr.append(arr[j,:].reshape(1,ROW_SIZE))
        i+=1
    f.close()
if __name__ == '__main__':
    #inputs
    coordinates = sys.argv[1] if len(sys.argv)>1 else 'coordinates.npy'
    #outputs
    combined = sys.argv[2] if len(sys.argv)>2 else 'coordinates.dbm'
    combine_coords(coordinates,combined)
