# Script to find the short axis distance in a certain mode
import numpy as np
import matplotlib.pyplot as plt
import sys
from biopandas.pdb import PandasPdb
from modes import plot_involvement, save_matrix, get_movie_muls
from gmx_file_processing import read_ndx, match_col_in_int_list

# shift the original pdb by a certain amplitude
# calculate the short axis distance of this mode
# store the short axis distance
# plot the short axis distance against modes to figure out peaks

# atom numbers of short axis
ALA63 = 491
PHE28 = 212

# atom indices in eigenmatrix
ALA63_i = ALA63 - 1
PHE28_i = PHE28 - 1

# get the short axis distance
def get_short_axis_distance(atom1, atom2):
    x2 = np.square(atom1[0] - atom2[0])
    y2 = np.square(atom1[1] - atom2[1])
    z2 = np.square(atom1[2] - atom2[2])
    return np.sqrt(x2 + y2 + z2)

# find atom coordinates based on atom number
def get_atom_coord(df, atom_number):
    stat_items = ['x_coord', 'y_coord', 'z_coord']
    coord = df[df['atom_number'] == atom_number] \
            [stat_items].to_numpy()[0]
    return coord

# calculate short axis distance of one mode
def get_one_mode_short_axis_delta_dist(short_axis_atoms_init, mode, indices, mul):
    # find atom mode vector by index (1X3)
    ALA_vector = mode[indices[0], :]
    PHE_vector = mode[indices[1], :]
    # harmonic shifting
    sample_steps = 150
    muls = get_movie_muls(mul, sample_steps)
    
    all_dist = []
    # find all distances
    for mul in muls: 
        # shift two atoms by +- mul
        ALA_coords = short_axis_atoms_init[0] + mul * ALA_vector
        PHE_coords = short_axis_atoms_init[1] + mul * PHE_vector
        all_dist.append(get_short_axis_distance(ALA_coords, PHE_coords))
    
    all_dist = np.array(all_dist)
    # find Dmax - Dmin
    delta_D = all_dist.max() - all_dist.min()
    
    return delta_D

# get short axis distances of all modes
def get_all_modes_short_axis_delta_dist(short_axis_atoms_init, eigenmatrix, eigenvalues, indices):
    # reshape eigenmatrix
    shift_shape = (int(eigenmatrix.shape[1]/3), 3)
    # list to store all short axis distances
    short_axis_delta_D = []
    for i in range(eigenmatrix.shape[0]):
        mode = eigenmatrix[:, i].reshape(shift_shape)
        mul = get_amp_from_ev(eigenvalues, i)
        one_delta_D = get_one_mode_short_axis_delta_dist(short_axis_atoms_init, mode, indices, mul)
        short_axis_delta_D.append(one_delta_D)
    
    return np.array(short_axis_delta_D)

# find difference in distances
def get_delta_D(D, D0):
    return D - D0
    
# plotting function
def plot_short_axis_spectrum(delta_D, title):
    plot_involvement(delta_D, title, mode_end=40, style='bars')

# get amplitude from eigenvalues
def get_amp_from_ev(ev, mode_index):
    eigenvalue = ev[mode_index] # A^2
    return np.sqrt(eigenvalue)

if __name__ == "__main__":
    # input pdb
    start_pdb = sys.argv[1]

    # read in as BioPandas Data Frame
    ppdb_start = PandasPdb()
    ppdb_start.read_pdb(start_pdb)

    # eigenmatrix
    eigenmatrix = np.load(sys.argv[2])
    
    # eigenvalues
    eigenvalues = np.load(sys.argv[3])

    # output figure filename
    if len(sys.argv) > 4:
        output_graph = sys.argv[4]
    else:
        output_graph = "bottleneck size increment [Å]"

    # shift amplitude
    # shift_amp = 40

    # short axis atom list
    ALA_init = get_atom_coord(ppdb_start.df['ATOM'], ALA63)
    PHE_init = get_atom_coord(ppdb_start.df['ATOM'], PHE28)
    short_axis_atoms = np.array([ALA_init, PHE_init])

    # initial_short_axis_distance
    D0 = get_short_axis_distance(ALA_init, PHE_init)

    # short axis indices
    sa_indices = np.array([ALA63_i, PHE28_i])
    # find Dmax
    delta_D = get_all_modes_short_axis_delta_dist(short_axis_atoms, eigenmatrix, eigenvalues, sa_indices)

    # store delta D
    # save_matrix("delta_short_axis_distance.npy", delta_D)

    # plot the spectrum
    plot_short_axis_spectrum(delta_D, output_graph)
    