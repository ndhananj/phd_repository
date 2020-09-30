# Script to find the short axis distance in a certain mode
import numpy as np
import matplotlib.pyplot as plt
import sys
from biopandas.pdb import PandasPdb
from modes import shift_by_mode, match_col_in_int_list

# shift the original pdb by a fixed amplitude
# calculate the short axis distance of this mode
# store the short axis distance
# plot the short axis distance against modes to figure out peaks

# get the short axis distance
def get_short_axis_distance(atom1, atom2):
    x2 = np.square(atom1[0] - atom2[0])
    y2 = np.square(atom1[1] - atom2[1])
    z2 = np.square(atom1[2] - atom2[2])
    return np.sqrt(x2 + y2 + z2)

# plotting function
def plot_short_axis_spectrum(delta_D, title):
    pass

# find atom coordinates based on atom number
def get_atom_coord(df, atom_number):
    stat_items = ['x_coord', 'y_coord', 'z_coord']
    coord = df['ATOM'][df['ATOM']['atom_number'] == atom_number] \
            [stat_items].to_numpy()[0]
    return coord

if __name__ == "__main__":
    # input pdb
    start_pdb = sys.argv[1]

    # read in as BioPandas Data Frame
    ppdb_start = PandasPdb()
    ppdb_start.read_pdb(start_pdb)

    # eigenmatrix
    eigenmatrix = sys.argv[2]

    # output figure filename
    if len(sys.argv) > 3:
        output_graph = sys.argv[3]
    else:
        output_graph = "short_axis_distance_across_all_modes"

    # atom numbers of short axis
    ALA63 = 491
    PHE28 = 212

    # shift amplitude
    shift_amp = 40

    ALA_init = get_atom_coord(ppdb_start.df, ALA63)
    PHE_init = get_atom_coord(ppdb_start.df, PHE28)

    # initial_short_axis_distance
    D0 = get_short_axis_distance(ALA_init, PHE_init)