from autocorrelation import *
from modes import *
import os.path
from os import path

def combine_means(means='means.npy',chunk_nums='chunk_nums.npy'):
    nums_list=load_matrix(chunk_nums)
    N=np.sum(nums_list)
    i=0
    means_list=[]
    while path.exists(chunk_name(means,i)):
        print('processing ',chunk_name(means,i))
        means_list.append(load_matrix(chunk_name(means,i))*nums_list[i])
        i+=1
    save_matrix(means,np.sum(means_list,axis=0)/N)
if __name__ == '__main__':
    #inputs
    means = sys.argv[1] if len(sys.argv)>1 else 'means.npy'
    chunk_nums = sys.argv[2] if len(sys.argv)>2 else 'chunk_nums.npy'
    combine_means(means,chunk_nums)
