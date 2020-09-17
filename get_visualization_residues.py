################################################################################
#  Get the visualization residues for alignment picture
#  Usage: python <this_file_name>
#  Example: python get_visualization_residues.py
################################################################################

from biopandas.pdb import PandasPdb
import pandas as pd
import numpy as np

def get_vis_res(input_file,prefix, chain_set):
    ppdb = PandasPdb()
    ppdb.read_pdb(input_file)
    df = ppdb.df['ATOM']

    chain = df['chain_id'].to_numpy()

    vis_res_list = []
    for chain_key in chain_set.keys():
        res_set = chain_set[chain_key]
        chains = df.loc[chain==chain_key]
        resn = chains['residue_number'].to_numpy()
        vis_res = chains.loc[np.isin(resn,res_set)]
        print(vis_res[['residue_name', 'residue_number']])
        vis_res_list.append(vis_res)

    print(vis_res_list)
    vis_res_df = pd.concat(vis_res_list)
    print(vis_res_df)
    vis_res_pdb = PandasPdb()
    vis_res_pdb.df['ATOM'] = vis_res_df
    vis_res_pdb.to_pdb(path=prefix+'_vis_res.pdb', \
        records=['ATOM'], gz=False, append_newline=True)
    return vis_res_df, vis_res_pdb


chain_set_4hea  = {'4':[87, 38], 'H':[216, 249]}
get_vis_res('4hea_no_dummies.pdb', '4hea', chain_set_4hea)
chain_set_6rfr  = {'C':[144, 95], '1':[232, 297]}
get_vis_res('6rfr_realigned.pdb', '6rfr', chain_set_6rfr)
