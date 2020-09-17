
from autocorrelation import *
from modes import *

def get_autocorrelation(\
    coordinates='coordinates.npy',\
    eigenmatrix='eigenmatrix.npy', \
    start_mode=0, end_mode=9):
    coords = load_matrix(coordinates)
    S = load_matrix(eigenmatrix)
    for i in range(start_mode,end_mode+1):
        Qcorr=autocorrelate_transformed(coords,S[i])
        plot_autocorrelate(Qcorr,i)

if __name__ == '__main__':
    coordinates = sys.arg[1] if len(sys.argv)>1 \
        else 'coordinates.npy'
    eigenmatrix = sys.argv[2] if len(sys.argv)>2 \
       else 'eigenmatrix.npy'
    start_mode = sys.argv[3] if len(sys.argv)>3 \
       else 0
    end_mode = sys.argv[4] if len(sys.argv)>4 \
       else 9
    get_autocorrelation(coordinates,eigenmatrix)
