################################################################################
#  Get atom entrance helices after input of helices
#  Usage: python <this_file_name>
#  Example: python get_atom_entrance.py
################################################################################

from biopandas.pdb import PandasPdb
import pandas as pd
import numpy as np

# return the column and column as a list
def getColAsList(df,col):
    df_col = df.filter(items=[col])
    list_df_col = [x[0] for x in df_col.values]
    return df_col, list_df_col

# return a list of indeces
def getFilteredColIndeces(df,col,to_match):
    df_col, list_df_col = getColAsList(df,col)
    indeces = [i for i in range(len(list_df_col)) if list_df_col[i]==to_match]
    return df_col, list_df_col, indeces

# return a list of indeces
def getFilteredCols(df,col,to_match):
    df_col, list_df_col, indeces = getFilteredColIndeces(df,col,to_match)
    filtered = df.iloc[indeces]
    return df_col, list_df_col, indeces, filtered

# dereference df assumed to be a single number
def getNumAt(df,row,col):
    return int(df.loc[[row],[col]].to_numpy()[0][0])

# returns true iff numberical values are betweedn init and end inclusive
def dfColBetween(df,col, init, end):
    return np.logical_and(df[col].to_numpy()>=init,df[col].to_numpy()<=end)

# select where colums are betweeen
def selRowsWhereColBetween(df,col,init,end):
    return df.loc[lambda df : dfColBetween(df,col, init, end)]

ppdb = PandasPdb()
ppdb.read_pdb('6rfr.pdb')
record_names, list_record_names, helix_indeces, helices = \
    getFilteredCols(ppdb.df['OTHERS'],'record_name','HELIX')

helix_entries, list_helix_entries = getColAsList(helices,'entry')
entry_parts = [x.split() for x in list_helix_entries]
entry_columns = \
    ['serNum', 'helixID', \
    'initResName', 'initChainID', 'initSeqNum', \
    'endResName', 'endChainID', 'endSeqNum', \
     'helixClass', 'length']
entry_df = pd.DataFrame(entry_parts,columns=entry_columns)

chains, list_of_chains, chain_H_indeces, chain_H_helices = \
    getFilteredCols(entry_df,'initChainID','1')

atoms, list_of_atoms, chain_H_atom_indeces, chain_H_atoms = \
    getFilteredCols(ppdb.df['ATOM'], 'chain_id', '1')

print(chain_H_helices)
for helix_idx in chain_H_indeces:
    initResNum = getNumAt(chain_H_helices,helix_idx,'initSeqNum')
    endResNum = getNumAt(chain_H_helices,helix_idx,'endSeqNum')
    print(initResNum,endResNum)
    helix_atoms = selRowsWhereColBetween(\
        chain_H_atoms,'residue_number',initResNum,endResNum)
    helix_pdb = PandasPdb()
    helix_pdb.df['ATOM'] = helix_atoms
    helix_pdb.to_pdb(path='helix'+str(helix_idx+1)+'.pdb', records=['ATOM'], \
        gz=False, append_newline=True)
