################################################################################
# This script will create a pdb of atoms
#   all the same element based on cavity statistics
# usage: python <this_file> <element> <stats_file> <out_file>
# example: python centers.py 'H' ref_cavs_stats.csv q_site_centers.pdb
################################################################################

import os, sys
from biopandas.pdb import PandasPdb
import pandas as pd
import numpy as np

def get_centers(stats_file):
    df = pd.read_csv(stats_file,header=None)
    coords = df[[10,11,12]]
    return coords

def output_centers(centers,element, out_file, xlabel=10, ylabel=11, zlabel=12):
    df_len = len(centers.index)
    df_data={ \
       'line_idx':np.arange(df_len), \
       'record_name':['ATOM']*df_len, \
       'atom_number':np.arange(1,df_len+1), \
       'blank_1':[' ']*df_len, \
       'atom_name':[element]*df_len, \
       'alt_loc':[' ']*df_len, \
       'residue_name':['DUM']*df_len, \
       'blank_2':[' ']*df_len, \
       'chain_id':['X']*df_len, \
       'residue_number':[1]*df_len, \
       'insertion':[' ']*df_len, \
       'blank_3':['   ']*df_len, \
       'x_coord':centers[[xlabel]].to_numpy().reshape((df_len,)), \
       'y_coord':centers[[ylabel]].to_numpy().reshape((df_len,)), \
       'z_coord':centers[[zlabel]].to_numpy().reshape((df_len,)), \
       'occupancy':[0.50]*df_len, \
       'b_factor':[35.88]*df_len, \
       'blank_4':['      ']*df_len, \
       'segment_id':['X1']*df_len, \
       'element_symbol':[element]*df_len, \
       'charge':[0.0]*df_len}
    df = pd.DataFrame(data=df_data)
    pdb = PandasPdb()
    pdb.df['ATOM'] = df
    pdb.to_pdb(path=out_file, records=['ATOM'], gz=False, append_newline=True)
    return df

if __name__ == "__main__":
    element = sys.argv[1]
    stats_file = sys.argv[2]
    out_file = sys.argv[3]

    centers = get_centers(stats_file)
    output_centers(centers,element, out_file)
