
from autocorrelation import *
from modes import *

def get_autocorrelation(\
    coordinates='coordinates.npy',\
    eigenmatrix='eigenmatrix.npy', \
    start_mode=0, end_mode=9):
    coords = load_matrix(coordinates)
    S = load_matrix(eigenmatrix)
    Q=transformed_coords(coords-coords.mean(axis=0),S)
    save_matrix('transformed_coords.npy', Q)
    N=min(40000,int(Q.shape[0]/5))
    for i in range(int(start_mode),int(end_mode)+1):
        Qcorr=autocorr1D(Q[:,i],N)
        save_matrix('autocorr_'+str(i)+'.npy',Qcorr[:N])
        plot_autocorrelate(Qcorr[:N],i)

if __name__ == '__main__':
    coordinates = sys.argv[1] if len(sys.argv)>1 \
        else 'coordinates.npy'
    eigenmatrix = sys.argv[2] if len(sys.argv)>2 \
       else 'eigenmatrix.npy'
    start_mode = sys.argv[3] if len(sys.argv)>3 \
       else 0
    end_mode = sys.argv[4] if len(sys.argv)>4 \
       else 9
    get_autocorrelation(coordinates,eigenmatrix,start_mode,end_mode)
