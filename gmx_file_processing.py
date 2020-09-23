################################################################################
# module to process files of various forms
# original author: Nithin Dhananjayan (ndhananj@ucdavis.edu)
################################################################################

import sys
import csv
import pandas as pd
import numpy as np

# This import registers the 3D projection, but is otherwise unused.
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401 unused import
import matplotlib.pyplot as plt

import re

import os

import glob

from biopandas.pdb import PandasPdb

from more_itertools import *

################################################################################
# General functions
################################################################################

# Return the first two group matches for the first search string that matches
#     It is assumed each search string produces two group matches
def return_2_match_groups(search_strings,string_to_search):
    f=lambda x:re.search(x,string_to_search)
    ms=[f(x) for x in search_strings]
    fulls = [ (m.group(1), m.group(2)) for m in ms if m ]
    if len(fulls[:])>0:
        return fulls[:][0][:]
    else:
        return [None, None]

# splice by list
def splice_by_idx_list(to_splice,splice_with):
    indeces = np.unique(splice_with,return_index=True)[1]
    idx_list = [splice_with[index] for index in sorted(indeces)]
    return [to_splice[i] for i in idx_list]

# take in lines that are sperated into sections and apply a function to it
def apply_func_lines_in_sections(f,secs):
    return [[f(line) for line in sec] for sec in secs]

# take in lines that are seperated into sections and split those lines
def split_lines_in_sections(secs):
    return apply_func_lines_in_sections(lambda line:line.split(),secs)

# given a set of data and columns, "harmonize" colums to the data
def harmonized_data_columns(data,cols):
    np_data = np.array(data)
    np_cols = np.array(cols).flatten()
    num_cols = np_data.shape[1]
    diff = num_cols-np_cols.shape[0]
    return np.pad(np_cols,diff)[diff:] if diff>0 else np_cols[:num_cols]

# wrapper to make data drame with data and parts using our defaults
def make_df_from_data_cols(data,cols):
    return pd.DataFrame(data=data,columns=harmonized_data_columns(data,cols))

# make a data frame from section columns and section section_parts
def make_df_from_sec_cols_parts(section_cols,section_parts):
    sections_range = range(len(section_cols))
    f = lambda i : make_df_from_data_cols(section_parts[i],section_cols[i])
    return [f(i) for i in sections_range]

# join the rows of a DataFrame
def join_rows_df(df):
    return df.astype(str).apply(lambda s: "    ".join(s),axis=1)

# join the rows of all DataFrames in a list
def join_rows_df_list(dfl):
    return [ join_rows_df(df) for df in dfl]

# join the rows of all DataFrames in a 2D list
def join_rows_df_list_2D(dfl):
    return [ join_rows_df_list(df) for df in dfl]

# concatenate lists
def concatenate_lists(l):
    return [item for sublist in l for item in sublist]

################################################################################
# data frame functions
################################################################################
def match_col_in_int_list(df,col,int_list):
    matches = np.array([df[col].to_numpy().astype(int) == i for i in int_list])
    return np.any(matches,axis=0)

def match_cols_in_int_list(df,cols,int_list):
    matches = np.array([match_col_in_int_list(df,col,int_list) for col in cols])
    return np.any(matches,axis=0)

def get_cols_in_int_list(df,cols,int_list):
    return df[match_cols_in_int_list(df,cols,int_list)]

def remap_integer_cols(df,cols,mapping):
    retVal = df.copy()
    for col in cols:
       idx_list = df[col].to_numpy().astype(int)
       retVal[col] = [mapping[i] for i in idx_list]
    return retVal


################################################################################
# pdb related functions
################################################################################
# the input is a biopandas read and a set of atom ids
def remove_pdb_atoms(pdb,atom_ids):
    remove_set = set(atom_ids)
    df = pdb.df['ATOM']
    col = df['atom_number'].to_numpy().astype(int)
    mark_f = np.vectorize(lambda x : x in remove_set)
    mark = np.invert(mark_f(col))
    pdb.df['ATOM'] = df[mark]
    return pdb

################################################################################
# xvg related functions
################################################################################
def process_xvg_params(params):
    search_strings = ['@(\w+)\s+(\w+)']
    search_strings.append('@\s+([\w\s]+)\s+\"(.+)\"')
    f=lambda x: return_2_match_groups(search_strings,x)
    param_pairs = [f(param) for param in params]
    y_axis_labels = [x[1] for x in param_pairs if x[0]=='yaxis  label']
    processed_params = {k:v for (k,v) in param_pairs }
    processed_params['yaxis  label'] = y_axis_labels
    return processed_params

def read_xvg_lines(filename):
    # read and split data into parts
    print("Reading xvg")
    return open(filename, 'r+').readlines()

def xvg_float_data_generator_from_lines(lines):
    return (np.array(line.split()).astype(np.float) \
        for line in lines if line[0] not in ['@','#'])

def get_xvg_data_array_from_lines(lines):
    return np.array(list(xvg_float_data_generator_from_lines(lines)))

def get_xvg_data_array_from_file(filename):
    return get_xvg_data_array_from_lines(open(filename, 'r+'))

def chunked_xvg_coord_data_from_file(filename, chunk_size=2000):
    return chunked(xvg_float_data_generator_from_lines(open(filename, 'r+')),\
        chunk_size)

def read_xvg(filename):
    lines = read_xvg_lines(filename)
    num_lines = len(lines)
    comments = (lines[i] for i in range(num_lines-1) if lines[i][0]=='#')
    params = (lines[i] for i in range(num_lines-1) if lines[i][0]=='@')
    print("xvg header processed")
    data_array = get_xvg_data_array_from_lines(lines,num_lines)
    # process and interpret the parts
    num_y_data_cols = data_array.shape[1]-1
    processed_params = process_xvg_params(params)
    y_labels = processed_params['yaxis  label']
    # use labels to place the data into a data frame
    num_y_labels = len(y_labels)
    num_new_labels = num_y_data_cols-num_y_labels
    if(num_new_labels>0): # if there aren't enough labels
        #print("Extending labels to "+str(num_y_data_cols)+" ...")
        new_labels = [y_labels[-1]+"_"+str(i) for i in range(num_new_labels)]
        y_labels.extend(new_labels)
        #print(y_labels)
    num_y_labels = len(y_labels)
    print("dependent variables found in xvg")
    # create dictionary to make data frame
    # start with y labels
    data_pairs = [(y_labels[i], data_array[:,i+1]) for i in range(num_y_labels)]
    # add the x label
    x_label = processed_params['xaxis  label']
    data_pairs.insert(0,(x_label,data_array[:,0]))
    # turn intoa dictionary
    data_dict = {k:v for (k,v) in data_pairs}
    print("data dictionary prepared")
    # turn into a data frame
    df = pd.DataFrame(data=data_dict)
    print("xvg interpreted")
    return {'data':df, 'xaxis label':x_label, 'yaxis labels':y_labels}

########
# assuming an xvg with an x_label and the rest as num_y_labels
# The input is the return value from reading the xvg

# returns  the average y_label values
def xvg_ylabel_avg(xvg_data):
    return xvg_data['data'][xvg_data['yaxis labels']].mean(axis=0)

# returns the first y_label values
def xvg_ylabel_first(xvg_data):
    return xvg_data['data'][xvg_data['yaxis labels']].iloc[0]

# returns the avg shift from first y_label values
def xvg_ylabel_shift(xvg_data):
    return xvg_ylabel_avg(xvg_data) - xvg_ylabel_first(xvg_data)

################################################################################
# Trajectory conversion functions
################################################################################

def convert_gro_to_pdb(file_prefix,regen=False):
    gro_file = file_prefix+'.gro'
    pdb_file = file_prefix+'.pdb'
    if(regen or not(os.path.exists(pdb_file))):
        cmd = \
          "echo '0' | gmx trjconv -f "+gro_file+' -s '+gro_file+' -o '+pdb_file
        os.system(cmd)

def get_frames_from_trj(trj_file, struct_file, beg, end, step, \
    out_prefix):
    cmd = "echo '0' | gmx trjconv -f "+trj_file+" -s "+struct_file+" -b "+\
        str(beg)+ " -e "+str(end)+" -dt "+\
        str(step)+" -sep -o "+out_prefix+".pdb"
    os.system(cmd)

def extract_position_from_traj_using_index(trj_file,struct_file,ndx_file,ndx,\
    output_file, beg=None, end=None, dt=None, regen=False):
    if(regen or not(os.path.exists(output_file))):
        cmd = "echo "+str(ndx)+" | gmx trajectory -f "+trj_file+" -s "+struct_file+\
            ' -n '+ndx_file+' -ox '+output_file
        if (beg is not None) and (end is not None) and (dt is not None) :
            cmd+=' -b '+beg+' -e '+end+' -dt '+dt
        os.system(cmd)

def extract_com_from_traj_using_index(trj_file,struct_file,ndx_file,ndx,\
    output_file,regen=False):
    if(regen or not(os.path.exists(output_file))):
        cmd = "echo "+str(ndx)+" | gmx traj -f "+trj_file+" -s "+struct_file+\
            ' -com -n '+ndx_file+' -ox '+output_file
        os.system(cmd)

def dump_snapshot(trj_file,struct_file,ndx_file,ndx,time,\
    output_file,regen=False):
    if(regen or not(os.path.exists(output_file))):
        cmd = "echo "+str(ndx)+" | gmx trjconv -f "+trj_file+\
            ' -dump '+str(time)+\
            " -s "+struct_file+ ' -n '+ndx_file+' -o '+output_file
        os.system(cmd)

################################################################################
# Simulation functions
################################################################################

#These functions assume the forcefield is set up to be available in directories

def grompp(mdp_file,init_struct_file,top_file,ndx_file,output_file):
    cmd = 'gmx grompp -f '+mdp_file+\
        ' -c '+init_struct_file+' -r '+init_struct_file+ \
        ' -p '+top_file+' -n '+ndx_file+' -o '+output_file
    os.system(cmd)

def mdrun(tpr_file,prefix,regen=False):
    if(regen or not(os.path.exists(prefix+'.xtc'))):
        cmd = 'gmx mdrun -s '+tpr_file+' -deffnm '+prefix
        os.system(cmd)

################################################################################
# Forcefield and Topology related functions
################################################################################

# make use of pdb2gmx create topologies and gro files
def make_neutral_peptide(input,choices):
    printf_str='\n'.join([str(c) for c in choices])
    [base,ext] = input.split('.')
    output = base+'.gro'
    top = base+'.top'
    itp = base+'_posre.itp'
    cmd = 'printf "' + printf_str + '" | gmx pdb2gmx -f '+ input + \
        ' -o ' + output + ' -p ' + top + ' -i ' + itp + \
        ' -chainsep id_or_ter -asp -ignh -posrefc 25000 -water none' + \
        ' -ff charmm36 -ter'
    print(cmd)
    os.system(cmd)

# convert a top to an itp by stripping forcfield and system information
def convert_top_to_itp(input):
    [base,ext]=input.split('.')
    output = base + '.itp'
    f = open(input, "r")
    lines=f.readlines()
    f.close()
    start_re = re.compile('\s*\[\s*moleculetype\s*\]')
    end_re = re.compile('\s*\[\s*system\s*\]')
    name_re = re.compile('\s*\;\s*Name')
    line_iter = range(len(lines))
    start = next(i for i in line_iter if start_re.match(lines[i]))
    end = next(i for i in line_iter if end_re.match(lines[i]))
    name_idx = next(i for i in line_iter if name_re.match(lines[i])) + 1
    lines[name_idx] = base+"            3\n"
    f2 = open(output,"w")
    f2.writelines(lines[start:end])
    f2.close()

# assumes there is at least 1 semicolon
def interpret_itp_comment_parts(parts):
    semicol_idx = [i for i in range(len(parts)) if parts[i]==';']
    semicol_idx.append(None)
    return parts[semicol_idx[0]+1:semicol_idx[1]]

# create a dictionary from headers and sections
def link_itp_sections_headers(headers,header_range,sections):
    section_for_header = {headers[j]:[] for j in header_range}
    for j in header_range:
        section_for_header[headers[j]].append(sections[j])
    return section_for_header

def read_itp(filename):
    # read and split data into parts
    f = open(filename, 'r+')
    # get lines
    lines = f.read().split("\n")
    num_lines = len(lines)
    lines_range = range(num_lines)
    # find the headers
    header_re = re.compile('\s*\[\s*(.*?)\s*\]')
    header_lines = [i for i in lines_range if header_re.match(lines[i])]
    headers = [header_re.match(lines[i]).group(1) for i in header_lines]
    num_headers = len(headers)
    header_range = range(num_headers)
    # Find the #ifdef or #endif section
    ifdef_re = re.compile('\#ifdef\s+([\w\_]+)')
    endif_re = re.compile('\#endif')
    ifdef_lines = [i for i in lines_range if ifdef_re.match(lines[i])]
    endif_lines = [i for i in lines_range if endif_re.match(lines[i])]
    ifdef_line_data = splice_by_idx_list(lines,ifdef_lines)
    endif_line_data = splice_by_idx_list(lines,endif_lines)
    ifdef_df = pd.DataFrame({'line_number':ifdef_lines,'line':ifdef_line_data})
    endif_df = pd.DataFrame({'line_number':endif_lines,'line':endif_line_data})
    # use headers to find sections
    end_func = lambda j: header_lines[j+1] if j+1 <num_headers else -1
    start_func = lambda j: header_lines[j]+1
    ends = [end_func(j) for j in header_range]
    starts = [start_func(j) for j in header_range]
    sections_meta_df = pd.DataFrame(\
       data={'headers':headers,'header_starts':starts,'header_ends':ends})
    # process sections
    raw_sections = [lines[starts[j]:ends[j]] for j in header_range]
    comment_line_re = re.compile('^\s*\;')
    blank_re = re.compile('^\s*$')
    to_ignore = [blank_re,ifdef_re,endif_re, comment_line_re]
    should_ignore = lambda line : np.any([e.match(line) for e in to_ignore])
    filter_sec = lambda sec : [line for line in sec if not(should_ignore(line))]
    filtered_secs = [ filter_sec(sec) for sec in raw_sections]
    # get the comment lines for each section
    should_include = lambda line : comment_line_re.match(line)
    sec_comments = lambda sec : [line for line in sec if should_include(line)]
    section_comment_lines = [ sec_comments(sec) for sec in raw_sections]
    section_comments = split_lines_in_sections(section_comment_lines)
    comments_df = pd.DataFrame(data={'comments':section_comment_lines})
    # put section information into data frames
    f = interpret_itp_comment_parts
    sec_cols = apply_func_lines_in_sections(f,section_comments)
    section_parts = split_lines_in_sections(filtered_secs)
    sections = make_df_from_sec_cols_parts(sec_cols,section_parts)
    section_for_header = \
        link_itp_sections_headers(headers,header_range,sections)
    return {'meta':sections_meta_df, **section_for_header, \
        'comments':comments_df, 'ifdef_lines':ifdef_df, 'endif_lines':endif_df,\
        'num_lines':num_lines}

def write_itp(itp_dict,filename):
    # set up output data structure
    lines = pd.DataFrame(columns=range(itp_dict['num_lines']),index=[0])
    # get data
    headers = itp_dict['meta']['headers'].to_numpy()
    comment_list = itp_dict['comments']['comments'].to_list()
    comments = ["".join(c) for c in comment_list]
    ifdef_lines = itp_dict['ifdef_lines']['line'].to_numpy()
    endif_lines = itp_dict['endif_lines']['line'].to_numpy()
    spliced_secs = splice_by_idx_list(itp_dict,headers)
    nested_secs=join_rows_df_list_2D(spliced_secs)
    sections = concatenate_lists(nested_secs)
    #get indeces
    starts = itp_dict['meta']['header_starts'].to_numpy()
    sec_starts = starts+2
    ends = itp_dict['meta']['header_ends'].to_numpy()
    num_sections = ends.shape[0]
    ends[num_sections-1]=itp_dict['num_lines']
    sec_ends = \
        [sec_starts[i]+sections[i].shape[0]-1 for i in range(num_sections)]
    true_sec_ends = np.min(np.stack([sec_ends,ends]),axis=0)
    sec_limits = np.stack([sec_starts,true_sec_ends])
    sec_ranges = [range(sec_limits[0,i], sec_limits[1,i]+1) \
        for i in range(num_sections)]
    ifdef_lns = itp_dict['ifdef_lines']['line_number'].to_numpy()
    endif_lns = itp_dict['endif_lines']['line_number'].to_numpy()
    # assign data into output data structure
    lines[starts]  = ['[ '+h+' ]' for h in headers]
    lines[starts+1] = comments
    lines[ifdef_lns] = "    ".join(ifdef_lines)
    lines[endif_lns] = "    ".join(endif_lines)
    for i in range(num_sections):
        lines[[j for j in sec_ranges[i]]] = sections[i].to_numpy()
    to_write = lines.T.fillna(' ')
    to_write.to_csv(filename,header=False,index=False,quotechar=" ")
    return to_write

# retrieve itp cols
def retrieve_itp_col(itp_dict,header,idx,col):
    return itp_dict[header][idx][col].to_numpy().astype(int)

# retrieve atom relevant cols and other data
def retrive_atom_relevant_itp(itp_dict):
    relevant_cols={'atoms':['nr'], 'bonds':['ai','aj'], 'pairs':['ai','aj'], \
        'angles':['ai','aj','ak'], 'dihedrals':['ai','aj','ak','al'], \
        'position_restraints':['atom']}
    headers = itp_dict['meta']['headers'].to_numpy()
    relevant_headers = [h for h in headers if h in relevant_cols.keys()]
    cols = [ relevant_cols[h] for h in relevant_headers]
    nested_secs = splice_by_idx_list(itp_dict,relevant_headers)
    sections = concatenate_lists(nested_secs)
    retrieved_cols = [ [sections[i][col].to_numpy().astype(int) \
        for col in cols[i]] for i in range(len(sections))]
    return relevant_cols, headers, relevant_headers, cols, nested_secs, \
        sections, retrieved_cols

# replace relevant headers
def replace_relavent_sections(itp_dict,relevant_headers,new_secs):
    header_range = range(len(relevant_headers))
    section_for_header=\
        link_itp_sections_headers(relevant_headers,header_range,new_secs)
    for h in relevant_headers:
        itp_dict[h]=section_for_header[h]
    return itp_dict

# apply column_func on columns and agregation function to make sections
def change_itp_cols_in_sections(col_f,agg_f,retrieved_cols):
    changed_cols = [ [col_f(col) for col in sec] for sec in retrieved_cols]
    changed_sections = [ agg_f(sec) for sec in changed_cols]
    return changed_cols, changed_sections

# take original sections and changed sections to produce full sections to store
def get_itp_sections_to_store(sections, changed_sections, f):
    return [f(sections,changed_sections,i) for i in range(len(sections))]

# return section with col replaced
def get_itp_replaced_section(cols,sec,csec,i):
    csec_df = pd.DataFrame(data=csec[i],columns=cols[i])
    sec[i]['new_idx'] = csec_df.index
    sec[i]=sec[i].set_index('new_idx')
    sec[i].loc[:,cols[i]] = csec_df.loc[:,cols[i]].copy()
    return sec[i]

# remove atoms (identified by number) from processed data structure
def remove_itp_atoms(itp_dict, atoms_to_remove):
    relevant_cols, headers, relevant_headers, cols, nested_secs, \
        sections, retrieved_cols = retrive_atom_relevant_itp(itp_dict)
    marking_func = \
        lambda col : np.any(np.array([col==a for a in atoms_to_remove]),axis=0)
    agg_func = lambda sec : np.any(np.array([cols for cols in sec]),axis=0)
    marked_cols, marked_sections = \
        change_itp_cols_in_sections(marking_func,agg_func,retrieved_cols)
    keep_f = lambda sec, csec, i: sec[i][np.invert(csec[i])]
    kept_sections = get_itp_sections_to_store(sections, marked_sections, keep_f)
    itp_dict = \
        replace_relavent_sections(itp_dict,relevant_headers,kept_sections)
    # renumber after
    relevant_cols, headers, relevant_headers, cols, nested_secs, \
        sections, retrieved_cols = retrive_atom_relevant_itp(itp_dict)
    orig_atom_numbers = itp_dict['atoms'][0]['nr'].to_numpy().astype(int)
    new_atom_numbers = \
       {orig_atom_numbers[i]:i+1 for i in range(len(orig_atom_numbers))}
    renumber = lambda x: new_atom_numbers[x] \
        if x in new_atom_numbers.keys() else x
    col_fv = np.vectorize(renumber)
    col_f = lambda col: col_fv(col) if col.shape[0]>0 else np.array([])
    agg_f = lambda sec: np.stack(sec).T
    renumbered_cols, renumbered_secs = \
       change_itp_cols_in_sections(col_f,agg_f,retrieved_cols)
    renumber_f = lambda sec, csec, i: \
        get_itp_replaced_section(cols,sec,csec,i)
    renumbered_sections = \
        get_itp_sections_to_store(sections, renumbered_secs, renumber_f)
    itp_dict = \
        replace_relavent_sections(itp_dict,relevant_headers,renumbered_sections)
    return itp_dict

# Set up monomer
def set_up_monomer(filename):
    [base,ext] = filename.split('.')
    mon_filename = base+'_mon.itp'
    cmd = "grep -v '; residue' " + filename + \
        " | grep -v '#' | grep -v '; Include' " + \
        " | perl -ne 's/\;\s*qtot.*$//; print' >" + mon_filename
    os.system(cmd)

# Fix termini charge in monomer
def fix_monomer_ter(itp_dict):
    ret_val = itp_dict.copy()
    ret_df = ret_val['atoms'][0].copy()
    nh1 = ret_df['nr']=='1'
    hn  = ret_df['nr']=='2'
    ct1 = ret_df['nr']=='4'
    ret_df.loc[nh1,'type']='NH1'
    ret_df.loc[nh1,'charge']=-0.47
    ret_df.loc[hn,'atom']='HN'
    ret_df.loc[hn,'charge']=0.31
    ret_df.loc[ct1,'charge']=0.07
    ht1 = ret_df['atom']=='HT1'
    ht2 = ret_df['atom']=='HT2'
    ret_df.loc[ht1,'type']='HA'
    ret_df.loc[ht2,'type']='HA'
    ret_df.loc[ht1,'charge']=0.09
    ret_df.loc[ht2,'charge']=0.09
    ret_val['atoms'][0]=ret_df
    return remove_itp_atoms(ret_val.copy(),[3])

# add inclusion of posre
def add_posre(file):
    [base,ext] = file.split('.')
    f=open(file,"a+")
    f.writelines(["\n",
        '; Include Position restraint file'+"\n",
        '#ifdef POSRES'+"\n",
        '#include "'+base+'_posre.itp"'+"\n",
        '#endif'+"\n"])
    f.close()

################################################################################
# ndx functions
################################################################################

def process_ndx_lines(lines_to_process):
    return [int(i) for i in \
        concatenate_lists([line.split() for line in lines_to_process])]

def read_ndx(filename):
    # read and split data into parts
    f = open(filename, 'r+')
    # get lines
    lines = f.read().split("\n")
    num_lines = len(lines)
    lines_range = range(num_lines)
    # find the headers
    header_re = re.compile('\s*\[\s*(.*?)\s*\]')
    header_lines = [i for i in lines_range if header_re.match(lines[i])]
    headers = [header_re.match(lines[i]).group(1) for i in header_lines]
    num_headers = len(headers)
    header_range = range(num_headers)
    retVal = {header:[] for header in headers}
    int_re = re.compile('\d+')
    for i in header_range:
        lines_to_process = lines[header_lines[i]+1:header_lines[i+1]] \
            if(i+1<num_headers) else lines[header_lines[i]+1:]
        retVal[headers[i]] = process_ndx_lines(lines_to_process)
    return retVal
