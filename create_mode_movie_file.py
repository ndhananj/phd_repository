from modes import *

def create_mode_movie_file(pdbfile,mode_idx,mul,ndxfile,ndx_name='System',\
    movie_steps=150,eigenmatrix="eigenmatrix.npy"):
    S = load_matrix(eigenmatrix)
    shift_shape=(int(S.shape[0]/3),3)
    mode = S[:,int(mode_idx)].reshape(shift_shape)
    mode_pdb_file=pdbfile+"_mode"+str(mode_idx)
    create_mode_movie(mul,movie_steps,\
        ndxfile,pdbfile,mode,mode_pdb_file,ndx_name)


if __name__ == '__main__':
    pdbfile = sys.argv[1]
    mode_idx = sys.argv[2]
    mul = sys.argv[3]
    ndxfile = sys.argv[4]
    ndx_name = sys.argv[5] if len(sys.argv) > 5 else 'System'
    movie_steps = sys.argv[6] if len(sys.argv) > 6 else 150
    eigenmatrix = sys.argv[7] if len(sys.argv) > 7 else "eigenmatrix.npy"
    create_mode_movie_file(pdbfile,mode_idx,mul,ndxfile,ndx_name,\
        movie_steps,eigenmatrix)
