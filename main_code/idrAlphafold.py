#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 10 10:26:27 2022

@author: magoncal
"""

import pandas as pd
import numpy as np
import re, os
import statistics as stats
from itertools import groupby

import resources

#comp_path = '/home/magoncal/Documents/data/projects/poly_cook/'
#un_prot = 'UP000005640'
#comp_path_un = os.path.join(comp_path, un_prot)

STRUC2DS_PAT = ['G', 'I', 'B', 'T', 'S', 'C']
STRUC2DS_REP = ['H', 'H', 'E', 'O', 'O', 'O']
STRUC2DS_DICT = {'H':0, 'E':1, 'O':2,' ':3, '|':4}

STRUC2D_DICT = {'H':0,'G':0,'I':0,'B':1,'E':1,'T':2,'S':2,'C':2,' ':3, '|':4}
STRUC2D_GRP = ['alpha-helix', 'beta-sheet','coil','unfolded', 'not aligned']


def extract_af_details(ss_data, ss_keys, ss_seq, std_sz=1400, std_over=200):
    ''' As alphafold fasta is too big to extract (90GB now) I'm using the
    CIF file output to contruct the sequence and size information. '''

    # Separating uniprot and models
    names_seqs = np.array(np.char.split(ss_keys, sep="_").tolist())
    models_int = np.array(names_seqs[:, 1], dtype="int64")
    base_names, idx_un, idx_inv, cnt_models = np.unique(names_seqs[:,0], return_index=True, return_inverse=True, return_counts=True)
    # Creating group
    base_rep = np.arange(0, len(idx_inv))[idx_inv]
    base_rep = np.column_stack((base_rep, models_int))
    # Reordering in decrescent order to make sure
    idx_desc = np.lexsort((-base_rep[:, 1], base_rep[:, 0]))
    keys_desc = ss_keys[idx_desc]
    keys_desc = keys_desc[idx_un]

    # Getting the sequence size using alphafold parameters
    coord_last = dict()
    for name, cnt in zip(keys_desc, cnt_models):
        if (cnt==1):
            seq_sz = len(ss_seq[name])
            coord_last[name] = [seq_sz, seq_sz, ss_seq[name], ss_data[name]]
        else:
            last_model = len(ss_seq[name])-std_sz-std_over
            # 1400 is alphafold limit now (std_sz)
            # -1 is for removing the last model,
            # 200 is the overlap made now on alphafold (std_over)
            # 1200 discounts the overlap since the start of the last part (std_sz-std_over)
            # +1 is required because its the size, not position
            seq_sz = std_sz+((cnt-2)*std_over)+last_model+1
            coord_last[name] = [seq_sz, last_model, ss_seq[name], ss_data[name]]
    
    # Generating a dataframe with the important information
    alphafold_data = pd.DataFrame.from_dict(coord_last, orient='index', columns=["af_len", "af_lastcoords", "af_seq", "af_ss"])
    alphafold_data["af_name"] = alphafold_data.index
    alphafold_data[["seq_name", "af_model"]] = alphafold_data.af_name.str.split("_", n=1, expand=True)
    alphafold_data = alphafold_data.reset_index(drop=True)
    alphafold_data = alphafold_data.loc[alphafold_data["af_model"]=="1", :]
    alphafold_data["af_ss_simpl"] = alphafold_data.af_ss.replace(STRUC2DS_PAT, STRUC2DS_REP, regex=True)
    alphafold_data["af_ss_simpl"] = alphafold_data.af_ss_simpl.replace("\s", "O")
    alphafold_data["af_model"] =  alphafold_data.af_model.astype("int64")
    # reordering columns
    cols = ['af_name', 'seq_name', 'af_model', 'af_len', 'af_lastcoords', 'af_seq', 'af_ss', 'af_ss_simpl']
    alphafold_data = alphafold_data[cols]
    return alphafold_data


def get_conf_scores(plddt_path, alphafold_data):
    ''' Gets the tab files with plddt scores and add it to a dataframe to be 
    sliced in the targetted regions. '''
    af_plddt = []
    with open(plddt_path, 'r') as handle:
        for line in enumerate(handle):
            rowsplit = line[1].rstrip("\n").split("\t")
            seq_name = rowsplit[0].split("_")[0]
            af_model = int(rowsplit[0].split("_")[1])
            scores = rowsplit[1].split(",")
            num_scores = [float(s) for s in scores]
            scores = ",".join(scores)
            af_plddt_coded = ""
            for v in num_scores:
                if (v < 50):
                    af_plddt_coded += "D"
                elif(v >= 50 and v < 70):
                    af_plddt_coded += "C"
                elif(v >= 70 and v < 90):
                    af_plddt_coded += "B"
                else:
                    af_plddt_coded += "A"
            lst_scores = [seq_name, af_model, scores, af_plddt_coded]
            af_plddt.append(lst_scores)
    cols = ["seq_name", "af_model", "af_plddt", "af_plddt_coded"]
    df_af_plddt = pd.DataFrame(af_plddt, columns=cols)
    alphafold_data = pd.merge(alphafold_data, df_af_plddt, how="left", on=["seq_name", "af_model"])
    return alphafold_data


def slice_scores(scores, start, end):
    ''' Gets just the region scores and calculate some statistics. '''
    num_scores = [float(s) for s in scores.split(",")][start:end]
    median_score = stats.median(num_scores)
    min_score = min(num_scores)
    reg_scores = ",".join([str(s) for s in num_scores])
    return pd.Series((v for v in [reg_scores, median_score, min_score]))


def filter_diff_seqs(ss_seq, seq_seq, alphafold_data):
    
    lst_remove = []
    
    for k, v in ss_seq.items():
        k_new = k.split("_")[0]
        if (k_new in seq_seq):
            if (v!=seq_seq[k_new]):
                lst_remove.append(k)
    
    alphafold_data = alphafold_data.loc[~alphafold_data.af_name.isin(lst_remove),:]
    return (alphafold_data)


def explode_surrounding_pLDDTs(polyXY_data_sel, delim, direct="left"):
    ''' Generate the matrix of 50 residues with pLDDT scores. '''
    
    scores = np.zeros((len(polyXY_data_sel), delim), dtype=float)
    
    # Get the values to a list
    list_scores = polyXY_data_sel.af_delim_left_plddt.str.split(",").tolist()
    # Get the size of each sublist and and generate a bool mask, reverting it 
    # for the left side of the sequence
    lens = [len(l) for l in list_scores]
    mask = np.arange(delim) < np.array(lens)[:,None]
    if (direct=="left"):
        mask = np.fliplr(mask)
    
    # Set the values to the right positions. In this cases np.concatenate flats
    # the list and the processing is faster.
    scores[mask] = np.concatenate(list_scores)
    return (scores)


def generate_surrounding_mask(polyXY_data_sel, delim, direct="left"):
    ''' Generate the matrix of 50 residues with poly and IDR mappings. '''
    pos_coords = np.zeros((len(polyXY_data_sel), delim), dtype=int)
    
    # 50 residue coordinates
    mask = np.arange(delim) < (np.array(polyXY_data_sel['af_delim_sz_'+direct].tolist()))[:,None]
    if (direct=="left"):
        mask = np.fliplr(mask)
    pos_coords[mask] = 1
    
    # IDR coordinates
    mask = np.arange(delim) < (np.array(polyXY_data_sel['af_delim_coord_idr_'+direct].tolist())+1)[:,None]
    if (direct=="left"):
        mask = np.fliplr(mask)
    pos_coords[mask] = 2
    
    # Poly coordinates
    mask = np.arange(delim) < (np.array(polyXY_data_sel['af_delim_coord_poly_'+direct].tolist())+1)[:,None]
    if (direct=="left"):
        mask = np.fliplr(mask)
    pos_coords[mask] = 3
    
    return(pos_coords)


def generate_mappings(polyXY_data_sel, comp_path_un, un_prot, source, delim, target):
    ''' Run the coordinates functions and save to disk in separate files. '''
    scores_left = explode_surrounding_pLDDTs(polyXY_data_sel, delim, "left")
    scores_right = explode_surrounding_pLDDTs(polyXY_data_sel, delim, "right")
    scores_100aa = np.hstack([scores_left, scores_right])
    if (target=="IDR"):
        target=""
    else:
        target="_"+target
    save_path = os.path.join(comp_path_un, un_prot+'_'+source+"_scores_"+str(delim)+"aa"+target+".csv")
    np.savetxt(save_path, scores_100aa, fmt="%03.2f", delimiter=",")
    
    map_left = generate_surrounding_mask(polyXY_data_sel, delim, "left")
    map_right = generate_surrounding_mask(polyXY_data_sel, delim, "right")
    coded_100aa = np.hstack([map_left, map_right])
    save_path = os.path.join(comp_path_un, un_prot+'_'+source+"_coded_"+str(delim)+"aa"+target+".csv")
    np.savetxt(save_path, coded_100aa, fmt="%d", delimiter=",")


def extract_ss_surround_child(polyXY_data, alphafold_data, delim=50):
    ''' Slices the surrounding region of the alphafold sequence and ss delim
    size to each side. As the procedure is the same, some extra fields will be
    added for slicing with IDR and poly.
    '''
    polyXY_data_sel = polyXY_data.loc[:, ["seq_name", "poly_name",  
                                          "poly_start", "poly_end", "poly_aa", "poly_size",
                                          "seq_len", "idr_start", "idr_end"]]
    polyXY_data_sel = pd.merge(polyXY_data_sel, alphafold_data, how="left", on="seq_name")
    polyXY_data_sel = polyXY_data_sel.loc[~polyXY_data_sel["af_name"].isna()&(polyXY_data_sel["poly_end"]<=polyXY_data_sel["af_len"]), :]
    
    # General cuts (poly) 
    polyXY_data_sel['poly_af_seq'] = polyXY_data_sel.apply(lambda x: x['af_seq'][int(x['poly_start']-1):int(x['poly_end'])], axis=1)
    polyXY_data_sel['poly_af_ss'] = polyXY_data_sel.apply(lambda x: x['af_ss'][int(x['poly_start']-1):int(x['poly_end'])], axis=1)
    polyXY_data_sel['poly_af_ss_simpl'] = polyXY_data_sel.apply(lambda x: x['af_ss_simpl'][int(x['poly_start']-1):int(x['poly_end'])], axis=1)
    polyXY_data_sel['poly_af_plddt_coded'] = polyXY_data_sel.apply(lambda x: x['af_plddt_coded'][int(x['poly_start']-1):int(x['poly_end'])], axis=1)
    polyXY_data_sel[['poly_af_plddt', 'poly_af_plddt_median', 'poly_af_plddt_min']] = polyXY_data_sel.apply(lambda x: slice_scores(x['af_plddt'], x['poly_start']-1, x['poly_end']), axis=1)
    
    # Extracting the delim surroundings
    polyXY_data_sel['af_ss_sz_left'] = (polyXY_data_sel['poly_size']*0.5).apply(np.floor).astype('int64')
    polyXY_data_sel['af_ss_sz_right'] = (polyXY_data_sel['poly_size']*0.5).apply(np.ceil).astype('int64')
    
    polyXY_data_sel['af_delim_center'] = polyXY_data_sel['poly_start']+polyXY_data_sel['af_ss_sz_left']
    polyXY_data_sel['af_delim_coord_poly_left'] = polyXY_data_sel['af_delim_center']-polyXY_data_sel['poly_start']
    polyXY_data_sel['af_delim_coord_poly_right'] = polyXY_data_sel['poly_end']-polyXY_data_sel['af_delim_center']+1
    polyXY_data_sel['af_delim_coord_idr_left'] = polyXY_data_sel['af_delim_center']-polyXY_data_sel['idr_start']-1
    polyXY_data_sel['af_delim_coord_idr_left'] = polyXY_data_sel.apply(lambda x: x['af_delim_coord_idr_left'] if x['af_delim_coord_idr_left']>=0 else 0, axis=1)
    polyXY_data_sel['af_delim_coord_idr_right'] = polyXY_data_sel['idr_end']-polyXY_data_sel['af_delim_center']+2
    polyXY_data_sel['af_delim_coord_idr_left'] = polyXY_data_sel.apply(lambda x: x['af_delim_coord_idr_right'] if x['af_delim_coord_idr_right']>=0 else 0, axis=1)
    
    polyXY_data_sel['af_align_left'] = polyXY_data_sel.apply(lambda x: x['af_seq'][:int(x['poly_start']-1+x['af_ss_sz_left'])], axis=1)
    polyXY_data_sel['af_align_right'] = polyXY_data_sel.apply(lambda x: x['af_seq'][int(x['poly_start']-1+x['af_ss_sz_left']):], axis=1)
    polyXY_data_sel['af_ss_align_left'] = polyXY_data_sel.apply(lambda x: x['af_ss'][:int(x['poly_start']-1+x['af_ss_sz_left'])], axis=1)
    polyXY_data_sel['af_ss_align_right'] = polyXY_data_sel.apply(lambda x: x['af_ss'][int(x['poly_start']-1+x['af_ss_sz_left']):], axis=1)
    polyXY_data_sel['af_ss_align_simpl_left'] = polyXY_data_sel.apply(lambda x: x['af_ss_simpl'][:int(x['poly_start']-1+x['af_ss_sz_left'])], axis=1)
    polyXY_data_sel['af_ss_align_simpl_right'] = polyXY_data_sel.apply(lambda x: x['af_ss_simpl'][int(x['poly_start']-1+x['af_ss_sz_left']):], axis=1)
    polyXY_data_sel['af_delim_sz_left'] = polyXY_data_sel.af_align_left.str.len()
    polyXY_data_sel['diff_left'] = delim-polyXY_data_sel.af_align_left.str.len()
    polyXY_data_sel['af_delim_sz_right'] = polyXY_data_sel.af_align_right.str.len()
    polyXY_data_sel['diff_right'] = delim-polyXY_data_sel.af_align_right.str.len()
    polyXY_data_sel['af_delim_left'] = polyXY_data_sel.apply(lambda x: '|'*x['diff_left']+x['af_align_left'] if x['diff_left']>0 else x['af_align_left'][abs(x['diff_left']):], axis=1)
    polyXY_data_sel['af_pos_delim_left_ori'] = polyXY_data_sel.poly_start+polyXY_data_sel.af_ss_sz_left-delim
    polyXY_data_sel['af_pos_delim_left'] = polyXY_data_sel.apply(lambda x: x['af_pos_delim_left_ori'] if x['af_pos_delim_left_ori']>0 else 1, axis=1)
    polyXY_data_sel['af_delim_right'] = polyXY_data_sel.apply(lambda x: x['af_align_right']+'|'*x['diff_right'] if x['diff_right']>0 else x['af_align_right'][:delim], axis=1)
    polyXY_data_sel['af_pos_delim_right'] = (polyXY_data_sel.poly_end+polyXY_data_sel.af_ss_sz_right)+delim
    polyXY_data_sel['af_pos_delim_right'] = polyXY_data_sel.apply(lambda x: x['af_pos_delim_right'] if x['af_pos_delim_right']<= x['seq_len'] else x['seq_len'], axis=1)
    polyXY_data_sel['af_ss_delim_left'] = polyXY_data_sel.apply(lambda x: '|'*x['diff_left']+x['af_ss_align_left'] if x['diff_left']>0 else x['af_ss_align_left'][abs(x['diff_left']):], axis=1)
    polyXY_data_sel['af_ss_delim_right'] = polyXY_data_sel.apply(lambda x: x['af_ss_align_right']+'|'*x['diff_right'] if x['diff_right']>0 else x['af_ss_align_right'][:delim], axis=1)
    polyXY_data_sel['af_ss_delim_simpl_left'] = polyXY_data_sel.apply(lambda x: '|'*x['diff_left']+x['af_ss_align_simpl_left'] if x['diff_left']>0 else x['af_ss_align_simpl_left'][abs(x['diff_left']):], axis=1)
    polyXY_data_sel['af_first_helix_left'] = polyXY_data_sel.apply(lambda x: first_ss(x['af_ss_delim_simpl_left'], "H", True), axis=1)
    polyXY_data_sel['af_first_sheet_left'] = polyXY_data_sel.apply(lambda x: first_ss(x['af_ss_delim_simpl_left'], "E", True), axis=1)
    #polyXY_data_sel['af_scr_helix_left', 'af_lastscr_helix_left'] = polyXY_data_sel.apply(lambda x: scr_region(x['af_ss_delim_simpl_left'], "H", True), axis=1)
    #polyXY_data_sel['af_scr_sheet_left', 'af_lastscr_sheet_left'] = polyXY_data_sel.apply(lambda x: scr_region(x['af_ss_delim_simpl_left'], "E", True), axis=1)
    polyXY_data_sel['af_ss_delim_simpl_right'] = polyXY_data_sel.apply(lambda x: x['af_ss_align_simpl_right']+'|'*x['diff_right'] if x['diff_right']>0 else x['af_ss_align_simpl_right'][:delim], axis=1)
    polyXY_data_sel['af_first_helix_right'] = polyXY_data_sel.apply(lambda x: first_ss(x['af_ss_delim_simpl_right'], "H", False), axis=1)
    polyXY_data_sel['af_first_sheet_right'] = polyXY_data_sel.apply(lambda x: first_ss(x['af_ss_delim_simpl_right'], "E", False), axis=1)
    #polyXY_data_sel['af_scr_helix_right', 'af_lastscr_helix_right'] = polyXY_data_sel.apply(lambda x: scr_region(x['af_ss_delim_simpl_right'], "H", False), axis=1)
    #polyXY_data_sel['af_scr_sheet_right', 'af_lastscr_sheet_right'] = polyXY_data_sel.apply(lambda x: scr_region(x['af_ss_delim_simpl_right'], "E", False), axis=1)
    
    # Had to calculate the scores first because of the transformations of the scores between str and int
    polyXY_data_sel['coords_left_start'] = polyXY_data_sel.apply(lambda x: abs(delim-(x['poly_start']-1+x['af_ss_sz_left'])) if delim-(x['poly_start']-1+x['af_ss_sz_left']) <0 else 0, axis=1)
    polyXY_data_sel['coords_left_end'] = polyXY_data_sel['poly_start']-1+polyXY_data_sel['af_ss_sz_left']
    polyXY_data_sel['coords_right_start'] = polyXY_data_sel['poly_start']-1+polyXY_data_sel['af_ss_sz_left']
    polyXY_data_sel['coords_right_end'] = polyXY_data_sel.apply(lambda x: x['poly_start']-1+x['af_ss_sz_left']+delim if x['poly_start']-1+x['af_ss_sz_left']+delim <= x['seq_len'] else x['seq_len'], axis=1)
    polyXY_data_sel['af_delim_left_plddt_coded'] = polyXY_data_sel.apply(lambda x: x['af_plddt_coded'][x['coords_left_start']:x['coords_left_end']], axis=1)
    polyXY_data_sel[['af_delim_left_plddt', 'af_delim_left_plddt_median', 'af_delim_left_plddt_min']] = polyXY_data_sel.apply(lambda x: slice_scores(x['af_plddt'], x['coords_left_start'], x['coords_left_end']), axis=1)
    polyXY_data_sel['af_delim_right_plddt_coded'] = polyXY_data_sel.apply(lambda x: x['af_plddt_coded'][x['coords_right_start']:x['coords_right_end']], axis=1)
    polyXY_data_sel[['af_delim_right_plddt', 'af_delim_right_plddt_median', 'af_delim_right_plddt_min']] = polyXY_data_sel.apply(lambda x: slice_scores(x['af_plddt'], x['coords_right_start'], x['coords_right_end']), axis=1)
    
    cols = ['poly_name', 'poly_size', 'poly_start', 'af_name', 
            'af_len','af_seq', 'af_ss', 'af_ss_simpl', 'af_plddt_coded', 
            'af_plddt', 'poly_af_seq', 'poly_af_ss', 'poly_af_ss_simpl', 
            'poly_af_plddt_coded', 'poly_af_plddt_median', 'poly_af_plddt_min', 
            'poly_af_plddt', 'af_delim_left',  'af_pos_delim_left_ori', 
            'af_pos_delim_left', 'af_ss_delim_left', 'af_ss_delim_simpl_left', 
            'af_first_helix_left', 'af_first_sheet_left', 'af_delim_left_plddt_coded', 
            'af_delim_left_plddt_median', 'af_delim_left_plddt_min', 'af_delim_left_plddt', 
            'af_delim_right', 'af_pos_delim_right', 'af_ss_delim_right', 'af_ss_delim_simpl_right', 
            'af_first_helix_right', 'af_first_sheet_right', 'af_delim_right_plddt_coded', 
            'af_delim_right_plddt_median', 'af_delim_right_plddt_min', 'af_delim_right_plddt',
            'af_delim_coord_poly_left', 'af_delim_coord_poly_right', 
            'af_delim_coord_idr_left', 'af_delim_coord_idr_right', 
            'af_delim_sz_left', 'af_delim_sz_right']
    polyXY_data_sel = polyXY_data_sel.loc[:, cols]
    return polyXY_data_sel


def extract_ss_surround(idr_data, alphafold_data):
    ''' Slices the surrounding region of the alphafold sequence and ss delim
    size to each side. As the procedure is the same, some extra fields will be
    added for slicing with IDR and poly.
    '''
    idr_data_sel = idr_data.loc[:, ["seq_name", "idr_name", "idr_start", "idr_end"]]
    idr_data_sel = idr_data_sel.astype({'idr_start': 'int32', 
                                              'idr_end': 'int32'})
    idr_data_sel = pd.merge(idr_data_sel, alphafold_data, how="left", on="seq_name")
    idr_data_sel = idr_data_sel.loc[~idr_data_sel["af_name"].isna(), :]
    
    # General cuts (IDR and poly)
    idr_data_sel['idr_af_seq'] = idr_data_sel.apply(lambda x: x['af_seq'][int(x['idr_start']-1):int(x['idr_end'])], axis=1)
    idr_data_sel['idr_af_ss'] = idr_data_sel.apply(lambda x: x['af_ss'][int(x['idr_start']-1):int(x['idr_end'])], axis=1)
    idr_data_sel['idr_af_ss_simpl'] = idr_data_sel.apply(lambda x: x['af_ss_simpl'][int(x['idr_start']-1):int(x['idr_end'])], axis=1)
    idr_data_sel['idr_af_plddt_coded'] = idr_data_sel.apply(lambda x: x['af_plddt_coded'][int(x['idr_start']-1):int(x['idr_end'])], axis=1)
    idr_data_sel[['idr_af_plddt', 'idr_af_plddt_median', 'idr_af_plddt_min']] = idr_data_sel.apply(lambda x: slice_scores(x['af_plddt'], x['idr_start']-1, x['idr_end']), axis=1)
    
    cols = ['idr_name', 'af_name', 'idr_af_seq', 'idr_af_ss', 'idr_af_ss_simpl', 
            'idr_af_plddt_coded', 'idr_af_plddt_median', 'idr_af_plddt_min', 
            'idr_af_plddt']
    idr_data_sel = idr_data_sel.loc[:, cols]
    return idr_data_sel


def count_dssp(dssp_pdbs, add_group=False):
    ''' Count the types of structures in the different provided region. '''
    struct=-1
    count_lst, prop_lst, map_lst, top_cont = [], [], [], []
    if (add_group):
        struct+=1
    for k, v in sorted(dssp_pdbs.items()):
        struc_lst = list(v)
        det_map = [STRUC2D_DICT[r] for r in struc_lst]
        det_count = np.bincount(det_map, minlength=len(STRUC2D_GRP)+struct)
        tot = np.sum(det_count)
        if tot>0:
            struc_prop = det_count/tot
        else:
            struc_prop = [0] * len(det_count)
        count_lst.append(det_count)
        if (add_group):
            map_lst.append(det_map)
        prop_lst.append(struc_prop)
    if (add_group):
        map_arr = np.array(map_lst)
    else:
        map_arr=np.empty([1,1])
    counts_arr = np.array(count_lst)
    prop_arr = np.array(prop_lst)
    return counts_arr, prop_arr, map_arr


def prepare_ss_counts(df_ss, prefix="idr", dir_name="", cut=""):
    ''' Extract the information of generic sets of columns with sequence and ss
    regions and format the output of the counting function. '''
    if (dir_name!=""):
        dir_name = "_"+dir_name
    if (cut!=""):
        cut = "_"+cut
    cols = list(df_ss.columns)
    # Dictionaries have no order, and zip creates a hash to order
    # To keep the order I'm forcing the order of the dictionary.
    idr_ss = dict(zip(df_ss[cols[0]], df_ss[cols[1]]))
    if (cut=='_delim'):
        ss_counts, ss_props, ss_map = count_dssp(idr_ss, True)
    else:
        ss_counts, ss_props, _  = count_dssp(idr_ss)
        ss_map=np.empty([1,1])
    
    # Merging the counts to add to the main dataframe
    ids_reg = df_ss[prefix+'_name'].values
    cols = [prefix+'_name', 'cnt_helix'+cut+dir_name, 
            'cnt_sheet'+cut+dir_name, 'cnt_coil'+cut+dir_name, 
            'cnt_unfolded'+cut+dir_name, 'prop_helix'+cut+dir_name, 
            'prop_sheet'+cut+dir_name, 'prop_coil'+cut+dir_name, 
            'prop_unfolded'+cut+dir_name]
    if (cut=='_delim'):
        cols.insert(5, 'cnt_noStruct'+cut+dir_name)
        cols.append('prop_noStruct'+cut+dir_name)
    df_dist2D = pd.DataFrame(np.hstack((ids_reg[:,None],
                                        ss_counts, ss_props)), columns=cols)
    df_dist2D = df_dist2D.astype({'cnt_helix'+cut+dir_name: 'int32', 
                                  'cnt_sheet'+cut+dir_name: 'int32',
                                  'cnt_coil'+cut+dir_name: 'int32',
                                  'cnt_unfolded'+cut+dir_name: 'int32',
                                  'prop_helix'+cut+dir_name: 'float32', 
                                  'prop_sheet'+cut+dir_name: 'float32',
                                  'prop_coil'+cut+dir_name: 'float32', 
                                  'prop_unfolded'+cut+dir_name: 'float32'})
    if (cut=='_delim'):
        df_dist2D = df_dist2D.astype({'cnt_noStruct'+cut+dir_name: 'int32', 
                                      'prop_noStruct'+cut+dir_name: 'float32'})
    df_dist2D['cnt_ss'+cut+dir_name] = df_dist2D.loc[:, ['cnt_helix'+cut+dir_name, 
                                                     'cnt_sheet'+cut+dir_name]].sum(axis=1)
    df_dist2D['prop_ss'+cut+dir_name] = df_dist2D.loc[:, ['prop_helix'+cut+dir_name, 
                                                      'prop_sheet'+cut+dir_name]].sum(axis=1)
    df_dist2D['cnt_other'+cut+dir_name] = df_dist2D.loc[:, ['cnt_coil'+cut+dir_name, 
                                                     'cnt_unfolded'+cut+dir_name]].sum(axis=1)
    df_dist2D['prop_other'+cut+dir_name] = df_dist2D.loc[:, ['prop_coil'+cut+dir_name, 
                                                     'prop_unfolded'+cut+dir_name]].sum(axis=1)
    return df_dist2D, ss_map


#simpl_ss = "O  O OO OOOOOOOOOOOO O         HHHHHHHHHHHHHHHHHHH"
#simpl_ss = "O  O OO OOOOOOOOOOOO O         OOOOOOOOOHHHHHHHHHH"
def scr_region(simpl_ss, struct, left):
    bin_arr = np.array([1 if l==struct else 0 for l in simpl_ss])
    if left:
        bin_arr = np.flip(bin_arr)
    cum_arr = np.cumsum(bin_arr)
    pos_arr = np.arange(1, len(cum_arr)+1)
    scr_arr = cum_arr/pos_arr
    return pd.Series((scr_arr, np.min(scr_arr)))


#k_mer=3
#ss = "HHOOO                                             "
def fuzzy_count_ss(k_mer, struct, ss, cnt_struct, left):
    ss_lst = np.array([1 if l==struct else 0 for l in ss])
    if left:
        ss_lst = np.flip(ss_lst)
    n = len(ss)
    xcol = np.arange(n-k_mer+1)
    inds_col = np.tile(np.arange(k_mer), (n-k_mer+1,1))
    inds_row = np.tile(xcol, (k_mer,1)).transpose()
    ss_int_slide = np.sum(ss_lst[inds_col + inds_row], axis=1)/k_mer
    cnt_struct = 1 if cnt_struct == 0 else cnt_struct #To avoid division by zero
    ss_max_idx = np.argmax(np.flip(ss_int_slide))/cnt_struct
    #ss_max_idx = len(ss_int_slide) if ss_max_idx==0 else ss_max_idx
    return ss_max_idx

#ss = "HOOOO                                             "
def first_ss(ss, struct, left):
    ss_lst = np.array([1 if l==struct else 0 for l in ss])
    if left:
        ss_lst = np.flip(ss_lst)
    ss_first_idx = np.argmax(ss_lst)+1
    ss_first_idx = 51 if all(ss_lst==0) else ss_first_idx
    #ss_max_idx = len(ss_int_slide) if ss_max_idx==0 else ss_max_idx
    return ss_first_idx

#test = df_idr_2D.loc[:, ['poly_name', 'af_ss_delim_simpl_left', 'af_ss_delim_left', 'cnt_helix_delim_left']]
#test['af_max3_helix_left'] = test.apply(lambda x: fuzzy_count_ss(3, "H", x['af_ss_delim_simpl_left'], x['cnt_helix_delim_left'], True), axis=1)

# Q93074_p_3_1
# af_pos_delim_left=1350
#ss = "O  OOOHHHHHHHHHHHHHHO OOOHHHHHHHHHHHHHHHHHHHHHHHHHOOOO        O   HHHHOOO OHHHHHHHHHOO HHHHHHHHHHHHH"
def struct_starts_ends(ss, af_pos_delim_left, struct):
    ''' Function extracts the coordinates of all helices or sheets from the provided
        region, returning the real sequence coords, size and number of units found 
        (turns or strands)
    '''
    ss_lst = np.array([1 if l==struct else 0 for l in ss])
    # Getting the indexes of all helices considering their position in the sequence
    idx_one = np.argwhere(ss_lst>0)+af_pos_delim_left
    # Getting the difference from each index. The regions with no idx are the separations between helices
    idx_diff = np.ediff1d(idx_one)!=1
    # Getting the start and end idx of each helix
    idx_helix = np.concatenate((idx_one[np.insert(idx_diff, 1, True)], idx_one[np.append(idx_diff, True)]), axis=1)
    # Getting the size of the helices
    idx_inter_diff = idx_helix[:,1]-idx_helix[:,0]
    # Ordering the indices to add them to the matrix
    idx_inter_ord = idx_inter_diff.argsort()
    idx_inter_diff = idx_inter_diff[idx_inter_ord, None]
    # Ordering the helices from smaller to bigger. This step is required to select always the biggest helix as tie breaker
    idx_helix = idx_helix[idx_inter_ord, :]
    # Adding the size and number of turns considering the canonical helical structure
    idx_helix = np.concatenate((idx_helix, idx_inter_diff), axis=1)
    return idx_helix.tolist()


def struct_select(poly_name, helices_coords, start, end, lim_sz):
    
    helices_coords = np.array(helices_coords)
    # Removing false helices. Sizes must be bigger than 3 initially, considering that helix3_10 are not found alone
    helices_coords = helices_coords[helices_coords[:,2]>lim_sz, :]
    # Getting the coordinates of the candidates to which the target structure belongs
    bool_pos = (end>=helices_coords[:, 0])&(start<=helices_coords[:, 1])
    # Get just the longest helix from the data
    helices_sel  = helices_coords[bool_pos, :]
    if (helices_sel.shape[0] > 0):
        helices_sel  = helices_coords[-1,:]
        return pd.Series(helices_sel.tolist())
    else:
        return pd.Series([np.nan, np.nan, np.nan])


def get_struct_details(df_2D_all, lim_sz, struct, struct_name, prefix):
    
    # Managing the Helix regions
    df_2D_struct = df_2D_all.loc[(df_2D_all['cnt_'+struct_name+'_delim_left']>lim_sz)|(df_2D_all['cnt_'+struct_name+'_delim_right']>lim_sz), :]
    df_2D_struct['af_ss_delim_simpl'] = df_2D_struct['af_ss_delim_simpl_left'].str.replace(" ", "-")+df_2D_struct['af_ss_delim_simpl_right'].str.replace(" ", "-")
    
    # Replaced the spaces by dashes to avoid errors
    df_2D_struct['af_'+struct_name+'_coords'] = df_2D_struct.apply(lambda x: struct_starts_ends(x['af_ss_delim_simpl'], x['af_pos_delim_left'], struct), axis=1)
    
    # Adding specific selection of best match of helix for the poly region
    cols = ['af_'+struct_name+'_'+prefix+'_start', 'af_'+struct_name+'_'+prefix+'_end', 'af_'+struct_name+'_'+prefix+'_size']
    df_2D_struct[cols] = df_2D_struct.apply(lambda x: struct_select(x['poly_name'], x['af_'+struct_name+'_coords'], x['poly_start'], x['poly_end'], lim_sz), axis=1)
    
    #teste['af_helix_poly_rel_start'] = teste['poly_start']/teste['af_helix_poly_start']
    #teste['af_helix_poly_rel_end'] = teste['poly_end']/teste['af_helix_poly_end']
    #teste['af_helix_idr_rel_start'] = teste['idr_start']/teste['af_helix_poly_start']
    #teste['af_helix_idr_rel_end'] = teste['idr_end']/teste['af_helix_poly_end']
    
    # Selecting just the target columns
    cols = ['poly_name', 'af_'+struct_name+'_coords']+cols
    df_2D_struct = df_2D_struct.loc[:, cols]
    
    return df_2D_struct   
    


def sel_valid(polyXY_data_sel, source, analysis_path):
    cols = ['poly_name', 'poly_size', 'poly_start', 'af_name', 
            'af_len', 'af_seq', 'af_ss', 'af_plddt_coded', 'poly_af_seq', 
            'poly_af_ss', 'poly_af_plddt_coded', 'af_delim_left', 'af_delim_right', 
            'af_ss_delim_left', 'af_ss_delim_right', 'af_delim_left_plddt_coded', 
            'af_delim_right_plddt_coded', 'af_plddt', 'poly_af_plddt', 
            'af_delim_left_plddt', 'af_delim_right_plddt']
    polyXY_data_valid = polyXY_data_sel.loc[:, cols]
    valid_name = 'EVAL_data_af_'+source+'.csv'
    polyXY_data_valid.to_csv(analysis_path+valid_name, index=False)
    
    
    
def final_cols(df_idr_2D):
    cols = ['poly_name', 'af_name', 'af_len', 'af_seq', 'af_ss', 'af_ss_simpl', 
            'af_plddt_coded', 'af_plddt', 'poly_af_seq', 'poly_af_ss', 'poly_af_ss_simpl', 
            'poly_af_plddt_coded', 'poly_af_plddt_median', 'poly_af_plddt_min', 
            'poly_af_plddt', 'cnt_helix_poly', 'cnt_sheet_poly', 'cnt_coil_poly', 
            'cnt_unfolded_poly', 'prop_helix_poly', 'prop_sheet_poly', 
            'prop_coil_poly', 'prop_unfolded_poly', 'cnt_ss_poly', 'prop_ss_poly', 
            'cnt_other_poly', 'prop_other_poly', 'af_delim_left', 'af_pos_delim_left_ori', 'af_pos_delim_left', 
            'af_ss_delim_left', 'af_ss_delim_simpl_left', 'af_first_helix_left', 
            'af_first_sheet_left', 'af_delim_left_plddt_coded', 'af_delim_left_plddt_median', 
            'af_delim_left_plddt_min', 'af_delim_left_plddt', 'cnt_helix_delim_left', 
            'cnt_sheet_delim_left', 'cnt_coil_delim_left', 'cnt_unfolded_delim_left', 
            'cnt_noStruct_delim_left', 'prop_helix_delim_left', 'prop_sheet_delim_left', 
            'prop_coil_delim_left', 'prop_unfolded_delim_left', 'prop_noStruct_delim_left', 
            'cnt_ss_delim_left', 'prop_ss_delim_left', 'cnt_other_delim_left', 
            'prop_other_delim_left', 'af_delim_right', 'af_pos_delim_right', 'af_ss_delim_right', 
            'af_ss_delim_simpl_right', 'af_first_helix_right', 'af_first_sheet_right', 
            'af_delim_right_plddt_coded', 'af_delim_right_plddt_median', 
            'af_delim_right_plddt_min', 'af_delim_right_plddt', 'cnt_helix_delim_right', 
            'cnt_sheet_delim_right', 'cnt_coil_delim_right', 'cnt_unfolded_delim_right', 
            'cnt_noStruct_delim_right', 'prop_helix_delim_right', 'prop_sheet_delim_right', 
            'prop_coil_delim_right', 'prop_unfolded_delim_right', 'prop_noStruct_delim_right', 
            'cnt_ss_delim_right', 'prop_ss_delim_right', 'cnt_other_delim_right', 
            'prop_other_delim_right', 'af_helix_coords', 'af_helix_poly_start',
            'af_helix_poly_end', 'af_helix_poly_size', 'af_sheet_coords',
            'af_sheet_poly_start', 'af_sheet_poly_end', 'af_sheet_poly_size']
    df_idr_2D = df_idr_2D.loc[:, cols]
    return(df_idr_2D)


def main_af(alphafold_ss_path, fasta_path, plddt_path, af_data_path):
    ss_data, ss_keys, ss_seq, _ = resources.extract_ss(alphafold_ss_path)
    _, _, seq_seq, _ = resources.extract_ss(fasta_path, s_type="seq", spl="|")
    alphafold_data = extract_af_details(ss_data, ss_keys, ss_seq)
    alphafold_data = get_conf_scores(plddt_path, alphafold_data)
    # Some of the Alphafold sequences are not in the same version as my Uniprots
    # so I need to remove them.
    alphafold_data = filter_diff_seqs(ss_seq, seq_seq, alphafold_data)
    alphafold_data.to_csv(af_data_path, index=False)
    

def main_idr(idr_path, af_data_path, comp_path_un, un_prot):
    alphafold_data = pd.read_csv(af_data_path)
    idr_data = pd.read_csv(idr_path, low_memory=False)
    idr_data = idr_data.sort_values(by=['idr_name'])
    idr_data_sel = extract_ss_surround(idr_data, alphafold_data)
    
    # idr
    df_ss = idr_data_sel.loc[:, ["idr_name", 'idr_af_ss', 'idr_af_seq']]
    df_idr_2D, _ = prepare_ss_counts(df_ss, prefix="idr", dir_name="", cut="idr")
    
    cols = ['idr_name', 'af_name', 'idr_af_seq', 'idr_af_ss', 
            'idr_af_ss_simpl', 'idr_af_plddt_coded', 'idr_af_plddt_median', 
            'idr_af_plddt_min', 'idr_af_plddt', 'cnt_helix_idr', 'cnt_sheet_idr',
            'cnt_coil_idr', 'cnt_unfolded_idr', 'prop_helix_idr', 'prop_sheet_idr',
            'prop_coil_idr', 'prop_unfolded_idr', 'cnt_ss_idr', 'prop_ss_idr', 
            'cnt_other_idr', 'prop_other_idr']

    df_idr_2D = pd.merge(idr_data_sel, df_idr_2D, how="left", on="idr_name", suffixes=('', '_y'))
    df_idr_2D = df_idr_2D.loc[:, cols]
    
    save_name = os.path.join(comp_path_un, un_prot+"idr_af_ss.csv")
    df_idr_2D.to_csv(save_name, index=False)
    

def main_poly(polyXY_path, af_data_path, comp_path_un, un_prot, source = "polyXY", valid = False, delim=50, lim_sz=3, target="IDR"):
    
    alphafold_data = pd.read_csv(af_data_path)
    polyXY_data = pd.read_csv(polyXY_path, low_memory=False)
    polyXY_data = polyXY_data.sort_values(by=['poly_name'])
    if target=="IDR":
        polyXY_data = polyXY_data.loc[~polyXY_data["idr_name"].isna(), :]
    else:
        polyXY_data = polyXY_data.loc[polyXY_data["idr_name"].isna(), :]
    
    # Duplication of poly_name can be caused by multiple IDR annotations. FIXED
    # Adding an extra ID based on poly_name to solve this problem.
    #polyXY_data["poly_idx"] = polyXY_data.groupby('poly_name').cumcount()+1
    #polyXY_data["poly_name"] = polyXY_data["poly_name"] + "_" + polyXY_data["poly_idx"].astype(str)
    
    polyXY_data_sel = extract_ss_surround_child(polyXY_data, alphafold_data)
    polyXY_data_sel = polyXY_data_sel.reset_index(drop=True)
    
    # Generate additional mappings
    generate_mappings(polyXY_data_sel, comp_path_un, un_prot, source, delim, target)
    
    # Make sure I don't rename my official file
    if (target=="IDR"):
        final_target = ""
        
    # poly
    df_ss = polyXY_data_sel.loc[:, ["poly_name", 'poly_af_ss', 'poly_af_seq']]
    df_idr_2D, _ = prepare_ss_counts(df_ss, prefix="poly", dir_name="", cut="poly")
    if (delim!=0):
        # left
        df_ss = polyXY_data_sel.loc[:, ["poly_name", 'af_ss_delim_left', 'af_delim_left']]
        df_cnt_af, ss_map_left = prepare_ss_counts(df_ss, prefix="poly", dir_name="left", cut="delim")
        df_idr_2D = pd.merge(df_idr_2D, df_cnt_af, how='left', on="poly_name")
        # right
        df_ss = polyXY_data_sel.loc[:, ["poly_name", 'af_ss_delim_right', 'af_delim_right']]
        df_cnt_af, ss_map_right = prepare_ss_counts(df_ss, prefix="poly", dir_name="right", cut="delim")
        df_idr_2D = pd.merge(df_idr_2D, df_cnt_af, how='left', on="poly_name")
    
    # Adding global sequence data in a group order (seq, poly, idr, left and right)
    df_idr_2D = pd.merge(polyXY_data_sel, df_idr_2D, how="left", on="poly_name", suffixes=('', '_y'))
    
    # Adding the structure annotations
    if (delim!=0):
        # Getting the extra IDR /poly info from the original dataframe
        cols = ['poly_name', 'poly_end']
        polyXY_data_sel2 = polyXY_data.loc[:, cols]
        
        df_2D_all = pd.merge(df_idr_2D, polyXY_data_sel2, how="left", on="poly_name")
        #df_2D_helix = df_2D_all.loc[]
        turn_sz = 3.6                           # https://www.sciencedirect.com/science/article/pii/0022283688906419?via%3Dihub
        strand_sz = 3 #??? (between 3 and 10)   # https://www.sciencedirect.com/science/article/pii/B9780128096338202677
        
        df_2D_helix = get_struct_details(df_2D_all, lim_sz, "H", "helix", "poly")
        df_2D_sheet = get_struct_details(df_2D_all, lim_sz, "E", "sheet", "poly")
    
        # Adding helix and extended specific data to dataframe
        df_idr_2D = pd.merge(df_idr_2D, df_2D_helix, how="left", on="poly_name")
        df_idr_2D = pd.merge(df_idr_2D, df_2D_sheet, how="left", on="poly_name")
    
        df_idr_2D = final_cols(df_idr_2D)
            
        # Mapping the SS by position
        map_delim_aa = np.hstack([ss_map_left, ss_map_right])
        save_path = os.path.join(comp_path_un, un_prot+'_'+source+"_af_"+str(delim)+"aa"+final_target+".csv")
        np.savetxt(save_path, map_delim_aa, fmt="%d", delimiter=",")
    
    if (valid):
        # The validation process is too slow with all columns. selecting just the ones that matter
        sel_valid(polyXY_data_sel, source, analysis_path)
    
    save_name = os.path.join(comp_path_un, un_prot+'_'+source+"_af_ss_"+str(delim)+"aa"+final_target+".csv")
    df_idr_2D.to_csv(save_name, index=False)


def count_top_contiguous(seq_aas):
    ''' Count the largest contiguous'''
    groups = sorted([list(g) for f, g in groupby(seq_aas)],key=len)
    last_group = groups[-1]
    before_group = groups[-2]
    if (len(last_group)==len(before_group)):
        res = "Equal"
    else:
        res = last_group[0]
    return pd.Series((res, len(last_group)))


def extra_polyXY(polyXY_path):
    df_poly_details = pd.read_csv(polyXY_path, low_memory=False)
    df_poly_details[['poly_long_aa','poly_long_cnt']] = df_poly_details.apply(lambda x: count_top_contiguous(x['poly_aa']), axis=1)
    df_poly_details.to_csv(polyXY_path, index=False)


def main_af_all(comp_path_un, un_prot):
    
    alphafold_ss_path = os.path.join(comp_path_un, un_prot+"_ss_alphafold.fasta")
    #polyXY_path = "/home/magoncal/Documents/data/projects/idr_cook/analysis/data_all_poly_coords.csv"
    polyXY_path = os.path.join(comp_path_un, un_prot+"_polyxy_details.csv")
    polyX_path = os.path.join(comp_path_un, un_prot+"_polyx_details.csv")
    idr_path = os.path.join(comp_path_un, un_prot+"_mobidb_lite_details.csv")
    plddt_path = os.path.join(comp_path_un, un_prot+"_plddts.tab")
    af_data_path = os.path.join(comp_path_un, un_prot+"_alphafold_data.csv")
    
    main_af(alphafold_ss_path, fasta_path, plddt_path, af_data_path)
    main_idr(idr_path, af_data_path, comp_path_un, un_prot)
    main_poly(polyX_path, af_data_path, comp_path_un, un_prot, "polyX", False, 50, 3)
    extra_polyXY(polyXY_path)
    main_poly(polyXY_path, af_data_path, comp_path_un, un_prot, "polyXY", False, 50, 3)
