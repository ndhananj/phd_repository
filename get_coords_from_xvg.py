from autocorrelation import *
from modes import *

def coords(xvg,coordfile):
    x = get_xvg_data_array_from_file(xvg)[:,1:]*10
    save_matrix(coordfile, x)

if __name__ == '__main__':
    xvg = sys.argv[1]
    coordfile = sys.argv[2]
    coords(xvg,coordfile)
