
from modes import *

import sys
import glob

def avg_cov_mat_stats(
    globPat='run*/covariance.npy', \
    covariance='avg/covariance.npy', \
    eigenvalues='avg/eigenvalues.npy', \
    eigenmatrix='avg/eigenmatrix.npy', \
    participations='avg/participations.npy'):
    filenames = [filename for filename in glob.glob(globPat)]
    cov = sum(load_matrix(filename) for filename in filenames)
    cov /= len(filenames)
    S, D, S2 = np.linalg.svd(cov, full_matrices=True)
    err = np.abs(S-S2)
    k=spring_constants_from_variances(D)
    P=get_all_atom_participations(S)
    print("cov ", cov.shape)
    print("D ", D.shape)
    print("S ", S.shape)
    print("S2 ",S2.shape)
    print("max|S-S2|=",np.max(err))
    print("mean(err)=",np.mean(err))
    save_matrix(covariance,cov)
    save_matrix(eigenvalues,D)
    save_matrix(spring_constants,k)
    save_matrix(eigenmatrix,S)
    save_matrix(participations,P)

if __name__ == '__main__':
    #inputs
    globPat = sys.argv[1] if len(sys.argv)>1 \
        else 'run*/covariance.npy'
    covariance = sys.arg[2] if len(sys.argv)>2 \
        else 'avg/covariance.npy'
    eigenvalues = sys.arg[3] if len(sys.argv)>3 \
        else 'avg/eigenvalues.npy'
    spring_constants = sys.argv[4] if len(sys.argv)> 4 \
       else 'avg/spring_constants.npy'
    eigenmatrix = sys.argv[5] if len(sys.argv)>5 \
       else 'avg/eigenmatrix.npy'
    participations = sys.argv[6] if len(sys.argv)>6 \
       else 'avg/participations.npy'
    # calculate
    avg_cov_mat_stats(
        globPat, \
        covariance, \
        eigenvalues, \
        eigenmatrix, \
        participations)
