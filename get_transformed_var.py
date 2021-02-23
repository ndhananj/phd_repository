
from autocorrelation import *
from modes import *
import os.path
from os import path

def get_transformed_var(\
    coordinates='transformed_coords.npy',\
    start_mode=0, end_mode=9,time_start=0
    ):
    coords = load_matrix(coordinates)
    for j in range(start_mode,end_mode+1):
        used = coords[:,j]
        variance = np.matmul(used.T,used)/used.shape[0]
        print("var(Q["+str(j)+"])=",variance)

if __name__ == '__main__':
    coordinates = sys.argv[1] if len(sys.argv)>1 \
        else 'transformed_coords.npy'
    start_mode = int(sys.argv[2]) if len(sys.argv)>2 \
       else 0
    end_mode = int(sys.argv[3]) if len(sys.argv)>3 \
       else 9
    time_start = float(sys.argv[4]) if len(sys.argv)>4 \
        else 0
    get_transformed_var(coordinates,start_mode,end_mode,time_start)
