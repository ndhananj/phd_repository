# Script to find the short axis distance in a certain mode
import numpy as np
import matplotlib.pyplot as plt
import sys
from biopandas.pdb import PandasPdb
from modes import shift_by_mode

if __name__ == "__main__":
    # input pdb
    start_pdb = sys.argv[1]

    # eigenmatrix
    eigenmatrix = sys.argv[2]

    # output figure filename
    if sys.argv[3]:
        output_graph = sys.argv[3]
    else:
        output_graph = "short_axis_distance_across_all_modes"