
from autocorrelation import *
from modes import *
import os.path
from os import path

def get_transformed_var(\
    coordinates='transformed_coords.npy',\
    start_mode=0, end_mode=9,time_start=0
    ):
    i=0
    names = []
    while path.exists(chunk_name(coordinates,i)):
        names.append(chunk_name(coordinates,i))
        i=i+1
    for j in range(start_mode,end_mode+1):
        coords = (load_matrix(name) for name in names)
        used = (c[:,j] for c in coords)
        var_pairs = ((np.matmul(u.T,u),u.shape[0]) for u in used)
        cum_var_pair = np.sum(np.array([[p[0],p[1]] for p in var_pairs]),axis=0)
        variance = cum_var_pair[0]/cum_var_pair[1]
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
