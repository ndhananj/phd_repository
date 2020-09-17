################################################################################
# This script will create a pdb for a box of atoms
#   all the same element spaced by spacing
# usage: python <this_file> <element> <spacing> \
#   <xmin> <xmax> <ymin> <ymax> <zmin> <zmax> <out_file>
# example: python box_of_element.py 'C' 1.7 5 40 -25 5 15 55 C_box.pdb
################################################################################

import os, sys
from biopandas.pdb import PandasPdb
import pandas as pd
import numpy as np

def box_of_element(element, spacing, \
    xmin, xmax, ymin, ymax, zmin, zmax, outfile):
    xs = np.arange(np.float(xmin),np.float(xmax),np.float(spacing))
    ys = np.arange(np.float(ymin),np.float(ymax),np.float(spacing))
    zs = np.arange(np.float(zmin),np.float(zmax),np.float(spacing))
    xx, yy, zz = np.meshgrid(xs,ys,zs)
    coords = np.stack([xx.flatten(), yy.flatten(), zz.flatten()], axis=1)
    len = coords.shape[0]
    df = pd.DataFrame(data={ \
       'line_idx':np.arange(len), \
       'record_name':['ATOM']*len, \
       'atom_number':np.arange(1,len+1), \
       'blank_1':[' ']*len, \
       'atom_name':[element]*len, \
       'alt_loc':[' ']*len, \
       'residue_name':['DUM']*len, \
       'blank_2':[' ']*len, \
       'chain_id':['X']*len, \
       'residue_number':[1]*len, \
       'insertion':[' ']*len, \
       'blank_3':['   ']*len, \
       'x_coord':coords[:,0], \
       'y_coord':coords[:,1], \
       'z_coord':coords[:,2], \
       'occupancy':[0.50]*len, \
       'b_factor':[35.88]*len, \
       'blank_4':['      ']*len, \
       'segment_id':['X1']*len, \
       'element_symbol':[element]*len, \
       'charge':[0.0]*len})
    pdb = PandasPdb()
    pdb.df['ATOM'] = df
    pdb.to_pdb(path=outfile, records=['ATOM'], gz=False, append_newline=True)
    return df

if __name__ == "__main__":
    element = sys.argv[1]
    spacing = sys.argv[2]

    xmin = sys.argv[3]
    xmax = sys.argv[4]
    ymin = sys.argv[5]
    ymax = sys.argv[6]
    zmin = sys.argv[7]
    zmax = sys.argv[8]

    outfile = sys.argv[9]

    df = box_of_element(element, spacing, \
        xmin, xmax, ymin, ymax, zmin, zmax, outfile)

    print(df)
