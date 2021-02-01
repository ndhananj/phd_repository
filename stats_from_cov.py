from modes import *

import sys

def get_cov_n_stats(
    covariance='covariance.npy', \
    eigenvalues='eigenvalues.npy', \
    eigenmatrix='eigenmatrix.npy', \
    participations='participations.npy'
    ):
    cov = load_matrix(covariance)
    S, D, S2 = svd(cov)
    err = np.abs(S-S2.T)
    k=spring_constants_from_variances(D)
    P=get_all_atom_participations(S)
    print("cov ", cov.shape)
    print("D ", D.shape)
    print("S ", S.shape)
    print("S2 ",S2.shape)
    print("max|S-S2.T|=",np.max(err))
    print("mean(err)=",np.mean(err))
    save_matrix(covariance,cov)
    save_matrix(eigenvalues,D)
    save_matrix(spring_constants,k)
    save_matrix(eigenmatrix,S)
    save_matrix(participations,P)

if __name__ == '__main__':
    #inputs
    covariance = sys.argv[1] if len(sys.argv)>1 \
        else 'covariance.npy'
        #outputs
    eigenvalues = sys.argv[2] if len(sys.argv)>2 \
        else 'eigenvalues.npy'
    spring_constants = sys.argv[3] if len(sys.argv)> 3 \
       else 'spring_constants.npy'
    eigenmatrix = sys.argv[4] if len(sys.argv)>4 \
       else 'eigenmatrix.npy'
    participations = sys.argv[5] if len(sys.argv)>5 \
       else 'participations.npy'
    # calculate
    get_cov_n_stats(
        covariance, \
        eigenvalues, \
        eigenmatrix, \
        participations
    )
