# Script to find the short axis distance in a certain mode
import numpy as np
import matplotlib.pyplot as plt
import sys
from biopandas.pdb import PandasPdb
from modes import shift_by_mode

# shift the original pdb by a fixed amplitude
# calculate the short axis distance of this mode
# store the short axis distance
# plot the short axis distance against modes to figure out peaks

# get the short axis distance
def get_short_axis_distance(atom1, atom2):
    pass

# plotting function
def plot_short_axis_spectrum():
    pass

# find atom coordinates based on atom number
def get_atom_coord(df, atom_number):
    coord = df['ATOM'][df['ATOM']['atom_number'] == atom_number]
    return coord

if __name__ == "__main__":
    # input pdb
    start_pdb = sys.argv[1]

    # read in as BioPandas Data Frame
    ppdb = PandasPdb()
    ppdb.read_pdb(start_pdb)

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

    ALA_coord = get_atom_coord(ppdb.df, ALA63)
    PHE_coord = get_atom_coord(ppdb.df, PHE28)