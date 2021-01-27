

from modes import *

import sys

def get_covariance_matrix_and_stats(xvg_input,\
    pdbToAlign=None, \
    outputForChunks=False,\
    unbias=False, \
    coordinates='coordinates.npy',\
    covariance='covariance.npy', \
    eigenvalues='eigenvalues.npy', \
    eigenmatrix='eigenmatrix.npy', \
    participations='participations.npy'):
    if(outputForChunks):
        resutGen = get_xvg_stats(xvg_input,fitfile=pdbToAlign,\
            outputForChunks=True,unbias=unbias)
        i=0
        for mean, cov, D, S, S2 in resutGen:
            err = np.abs(S-S2)
            k=spring_constants_from_variances(D)
            P=get_all_atom_participations(S)
            print("cov ", cov.shape)
            print("D ", D.shape)
            print("S ", S.shape)
            print("S2 ",S2.shape)
            print("max|S-S2|=",np.max(err))
            print("mean(err)=",np.mean(err))
            save_matrix(chunk_name(covariance,i),cov)
            save_matrix(chunk_name(eigenvalues,i),D)
            save_matrix(chunk_name(spring_constants,i),k)
            save_matrix(chunk_name(eigenmatrix,i),S)
            save_matrix(chunk_name(participations,i),P)
            i+=1
    else:
        mean, cov, D, S, S2, coords = \
            get_xvg_stats(xvg_input,fitfile=pdbToAlign,\
            outputForChunks=False,unbias=unbias)
        err = np.abs(S-S2)
        k=spring_constants_from_variances(D)
        P=get_all_atom_participations(S)
        print("coords ", coords.shape)
        print("cov ", cov.shape)
        print("D ", D.shape)
        print("S ", S.shape)
        print("S2 ",S2.shape)
        print("max|S-S2|=",np.max(err))
        print("mean(err)=",np.mean(err))
        save_matrix(coordinates,coords)
        save_matrix(covariance,cov)
        save_matrix(eigenvalues,D)
        save_matrix(spring_constants,k)
        save_matrix(eigenmatrix,S)
        save_matrix(participations,P)

if __name__ == '__main__':
    #inputs
    xvg_input = sys.argv[1]
    pdbToAlign = sys.argv[2] if len(sys.argv)>2 else None
    #parameters
    outputForChunks = bool(sys.argv[3]) if len(sys.argv)>3 else False
    unbias = bool(sys.argv[4]) if len(sys.argv)>4 else False
    #outputs
    coordinates = sys.argv[5] if len(sys.argv)>5 \
        else 'coordinates.npy'
    covariance = sys.argv[6] if len(sys.argv)>6 \
        else 'covariance.npy'
    eigenvalues = sys.argv[7] if len(sys.argv)>7 \
        else 'eigenvalues.npy'
    spring_constants = sys.argv[8] if len(sys.argv)> 8 \
       else 'spring_constants.npy'
    eigenmatrix = sys.argv[9] if len(sys.argv)>9 \
       else 'eigenmatrix.npy'
    participations = sys.argv[10] if len(sys.argv)>10 \
       else 'participations.npy'
    # calculate
    get_covariance_matrix_and_stats(xvg_input,\
        pdbToAlign, \
        outputForChunks, \
        unbias, \
        coordinates, \
        covariance, \
        eigenvalues, \
        eigenmatrix, \
        participations)
