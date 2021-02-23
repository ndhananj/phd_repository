from modes import *

import sys

def get_cov_n_stats(
    time_step=0.05, \
    time_start=1500,\
    coordinates='coordinates.npy',\
    covariance='covariance.npy', \
    eigenvalues='eigenvalues.npy', \
    eigenmatrix='eigenmatrix.npy', \
    participations='participations.npy'
    ):
    coords = load_matrix(coordinates)
    idx_start=int(time_start/time_step)
    print("Starting index is ", idx_start)
    print("Calculating stats...")
    stat_coords=coords[idx_start:,:]
    print("stat. coords",stat_coords.shape)
    mean, cov, D, S, S2 = calc_single_coord_stats(stat_coords)
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
    #parameters
    time_step = float(sys.argv[1]) if len(sys.argv)>1 else 0.05
    time_start = float(sys.argv[2]) if len(sys.argv)>2 else 1500
    #inputs
    coordinates = sys.argv[3] if len(sys.argv)>3 \
        else 'coordinates.npy'
    #outputs
    covariance = sys.argv[4] if len(sys.argv)>4 \
        else 'covariance.npy'
    eigenvalues = sys.argv[5] if len(sys.argv)>5 \
        else 'eigenvalues.npy'
    spring_constants = sys.argv[6] if len(sys.argv)> 6 \
       else 'spring_constants.npy'
    eigenmatrix = sys.argv[7] if len(sys.argv)>7 \
       else 'eigenmatrix.npy'
    participations = sys.argv[8] if len(sys.argv)>8 \
       else 'participations.npy'
    # calculate
    get_cov_n_stats(
        time_step,\
        time_start, \
        coordinates, \
        covariance, \
        eigenvalues, \
        eigenmatrix, \
        participations
    )
