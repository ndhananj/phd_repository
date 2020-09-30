# Script to find the short axis distance in a certain mode
import numpy as np
import matplotlib.pyplot as plt
import sys
from biopandas.pdb import PandasPdb
from modes import shift_by_mode
from gmx_file_processing import read_ndx

# shift the original pdb by a fixed amplitude
# calculate the short axis distance of this mode
# store the short axis distance
# plot the short axis distance against modes to figure out peaks

# atom numbers of short axis
ALA63 = 491
PHE28 = 212

# get the short axis distance
def get_short_axis_distance(atom1, atom2):
    x2 = np.square(atom1[0] - atom2[0])
    y2 = np.square(atom1[1] - atom2[1])
    z2 = np.square(atom1[2] - atom2[2])
    return np.sqrt(x2 + y2 + z2)

# find atom coordinates based on atom number
def get_atom_coord(df, atom_number):
    stat_items = ['x_coord', 'y_coord', 'z_coord']
    coord = df['ATOM'][df['ATOM']['atom_number'] == atom_number] \
            [stat_items].to_numpy()[0]
    return coord

# get short axis distances of all modes
def get_all_modes_short_axis_dist(df, eigenmatrix, indices, mul):
    short_axis_dist_list = []
    for i in range(eigenmatrix.shape[0]):
        shifted_df = shift_by_mode(df, eigenmatrix[:,i], indices, mul)
        ALA_i = get_atom_coord(shifted_df, ALA63)
        PHE_i = get_atom_coord(shifted_df, PHE28)
        D = get_short_axis_distance(ALA_i, PHE_i)
        short_axis_dist_list.append(D)
    return np.array(short_axis_dist_list)
    
# plotting function
def plot_short_axis_spectrum(delta_D, title):
    pass

if __name__ == "__main__":
    # input pdb
    start_pdb = sys.argv[1]

    # read in as BioPandas Data Frame
    ppdb_start = PandasPdb()
    ppdb_start.read_pdb(start_pdb)

    # index file
    ndx = read_ndx(sys.argv[2])

    # eigenmatrix
    eigenmatrix = sys.argv[3]

    # output figure filename
    if len(sys.argv) > 4:
        output_graph = sys.argv[4]
    else:
        output_graph = "short_axis_distance_across_all_modes"

    # shift amplitude
    shift_amp = 40

    ALA_init = get_atom_coord(ppdb_start.df, ALA63)
    PHE_init = get_atom_coord(ppdb_start.df, PHE28)

    # initial_short_axis_distance
    D0 = get_short_axis_distance(ALA_init, PHE_init)