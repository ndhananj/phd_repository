
from modes import *

def get_involvement(involved_residues,\
    outputForChunks=False,\
    participation='participations.npy',\
    residues='residues.npy',\
    involvement_string='Bottleneck Involvement'):
    if(outputForChunks):
        numChunks=len(chunk_glob(participation))
        for i in range(numChunks):
            get_involvement(involved_residues,False,\
                chunk_name(participation,i),residues,\
                chunk_name(involvement_string,i))
    else:
        P=load_matrix(participation)
        resi=load_matrix(residues)
        f = open(involved_residues, "r")
        toInclude = np.array(f.read().split()).astype(int)
        I=involvement_in_mode_based_on_participation(\
            P,resi,toInclude)
        save_matrix(involvement_string+'.npy',I)
        plot_involvement(I,involvement_string)

if __name__ == '__main__':
    involved_residues = sys.argv[1]
    outputForChunks = bool(sys.argv[2]) if len(sys.argv)>2 \
        else False
    participation = sys.argv[3] if len(sys.argv)>3 \
       else 'participations.npy'
    residues = sys.argv[4] if len(sys.argv)>4 \
       else 'residues.npy'
    involvement_string = sys.argv[5] if len(sys.argv)>5 \
       else 'Bottleneck Involvement'
    get_involvement(involved_residues,outputForChunks,participation,\
        residues,involvement_string)
