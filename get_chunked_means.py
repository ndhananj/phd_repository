from autocorrelation import *
from modes import *
import os.path
from os import path

def get_means(
    coordinates='coordinates.npy',means='means.npy',chunk_nums='chunk_nums.npy'
    ):
    i=0
    nums_list=[]
    while path.exists(chunk_name(coordinates,i)):
        print(chunk_name(coordinates,i))
        coords = load_matrix(chunk_name(coordinates,i))
        nums_list.append(coords.shape[0])
        save_matrix(chunk_name(means,i),coords.mean(axis=0))
        i+=1
    save_matrix(chunk_nums,np.array(nums_list))
if __name__ == '__main__':
    #inputs
    coordinates = sys.argv[1] if len(sys.argv)>1 else 'coordinates.npy'
    #outputs
    means = sys.argv[2] if len(sys.argv)>2 else 'means.npy'
    chunk_nums = sys.argv[3] if len(sys.argv)>3 else 'chunk_nums.npy'
    get_means(coordinates,means,chunk_nums)
