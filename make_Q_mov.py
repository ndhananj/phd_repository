
from autocorrelation import *
from modes import *

def make_Q_mov(\
    coordinates,ndxfile,pdbfile,newpdbfile,ndx_name,stride=1,\
    eignevector='eigenvector0.npy',time_start=0,time_step=0.05):
    coords = load_matrix(coordinates)
    print("coords",coords.shape)
    muls = coords[::stride]
    print("muls",muls.shape)
    v = load_matrix(eignevector)
    print("v",v.shape)
    natoms = int(v.shape[0]/3)
    mode = v.reshape((natoms,3))
    print("mode",mode.shape)
    make_movie_from_muls(muls,ndxfile,pdbfile,mode,newpdbfile,ndx_name)

if __name__ == '__main__':
    coordinates = sys.argv[1]
    ndxfile = sys.argv[2]
    pdbfile = sys.argv[3]
    newpdbfile = sys.argv[4]
    ndx_name = sys.argv[5]
    stride = int(sys.argv[6]) if len(sys.argv)>6 else 1
    mode = sys.argv[7] if len(sys.argv)>7 else 'eigenvector0.npy'
    time_start = float(sys.argv[8]) if len(sys.argv)>8 else 0
    time_step = float(sys.argv[9]) if len(sys.argv)>9 else 0.05
    make_Q_mov(\
        coordinates,ndxfile,pdbfile,newpdbfile,ndx_name,\
        stride,mode,time_start,time_step)
