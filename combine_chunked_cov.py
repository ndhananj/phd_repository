from autocorrelation import *
from modes import *
import os.path
from os import path

def combine_cov(
    covariance='covariance.npy',means='means.npy',chunk_nums='chunk_nums.npy'
    ):
    nums_list=load_matrix(chunk_nums)
    N=np.sum(nums_list)
    i=0
    addends=[]
    while \
    path.exists(chunk_name(covariance,i)) and path.exists(chunk_name(means,i)):
        print('processing ',chunk_name(means,i))
        mean = load_matrix(chunk_name(means,i))
        shifts = \
            mean.reshape(mean.shape[0],1).dot(mean.reshape(1,mean.shape[0]))
        print('processing ',chunk_name(covariance,i))
        cov = load_matrix(chunk_name(covariance,i))
        addends.append((cov+shifts)*nums_list[i])
        i+=1
    mean = load_matrix(means)
    shift = \
        mean.reshape(mean.shape[0],1).dot(mean.reshape(1,mean.shape[0]))
    save_matrix(covariance,np.sum(addends,axis=0)/N-shift)
if __name__ == '__main__':
    #inputs
    covariance = sys.argv[1] if len(sys.argv)>1 else 'covariance.npy'
    means = sys.argv[2] if len(sys.argv)>2 else 'means.npy'
    chunk_nums = sys.argv[3] if len(sys.argv)>3 else 'chunk_nums.npy'
    combine_cov(covariance,means,chunk_nums)
