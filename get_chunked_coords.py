from autocorrelation import *
from modes import *

def get_coords(\
    xvgfile,pdbToAlign=None,unbias=False,coordinates='coordinates.npy'
    ):
    coords = fitted_chunks_from_xvg(xvgfile,fitfile=pdbToAlign,unbias=unbias)
    i=0
    for c in coords:
        print(chunk_name(coordinates,i))
        save_matrix(chunk_name(coordinates,i),c)
        i+=1
if __name__ == '__main__':
    #inputs
    xvg_input = sys.argv[1]
    pdbToAlign = sys.argv[2] if len(sys.argv)>2 else None
    #parameters
    unbias = bool(sys.argv[3]) if len(sys.argv)>3 else False
    #outputs
    coordinates = sys.argv[5] if len(sys.argv)>5 \
        else 'coordinates.npy'
    get_coords(xvg_input,pdbToAlign, unbias, coordinates)
