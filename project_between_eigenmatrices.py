from modes import *

import sys

def project_between_eigenmatrices(dir1,dir2):
    e1=load_matrix(dir1+'/eigenmatrix.npy')
    e2=load_matrix(dir2+'/eigenmatrix.npy')
    proj=np.matmul(e2[:,:40].T,e1[:,:40])
    print(proj)
    fig, ax = plt.subplots()
    im = ax.imshow(proj)
    ax.set_xticks(np.arange(40))
    ax.set_yticks(np.arange(40))
    cbar = ax.figure.colorbar(im, ax=ax)
    plt.show()

if __name__ == '__main__':
    # first directory of eigenmatrix
    dir1 = sys.argv[1]
    # second directory of eigenmatrix
    dir2 = sys.argv[2]
    project_between_eigenmatrices(dir1,dir2)
