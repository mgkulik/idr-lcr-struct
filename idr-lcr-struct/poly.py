#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 15 15:23:01 2021

@author: magoncal
"""

from Bio import SeqIO, SearchIO
import os
import numpy as np
import pandas as pd
import pickle
import collections
import math

import resources
        

# NOT IN USE
def get_polyXY(polyXY_path, polytab_path, sep="\t"):
    ''' Transforms the original tab file with one polyXY per row in a list of
    coordinates with one protein por row. '''
    old_id = ""
    tab_rows = []
    uniprots = []
    with open(polyXY_path, 'r') as handle:
        next(handle)
        for line in enumerate(handle):
            rowsplit = line[1].rstrip("\n").split("\t")
            uniprot_id = rowsplit[0].split("|")[1]
            uniprots.append(uniprot_id)
            reg = str(rowsplit[1])+'-'+str(rowsplit[2])
            if old_id!=uniprot_id:
                if old_id!="":
                    row_content = old_id+sep+reg_str+'\n'
                    tab_rows.append(row_content)
                reg_str = reg
            else:
                reg_str = reg_str + sep + reg
            old_id = uniprot_id
    row_content = uniprot_id+sep+reg_str+'\n'
    tab_rows.append(row_content)
    resources.save_file(tab_rows, polytab_path)
    return uniprots
  

def extract_polyXY_features(polyXY_path, fasta_data, tp=2):
    ''' Extract the repeats from the tab file and add some additional properties. '''
    poly_lst=[]
    old_id=""
    with open(polyXY_path, 'r') as handle:
        next(handle)
        for line in enumerate(handle):
            rowsplit = line[1].rstrip("\n").split("\t")
            uniprot_id = rowsplit[0].split("|")[1]
            seq_aa = fasta_data[uniprot_id][0]
            seq_len = len(seq_aa)
            seq_desc = fasta_data[uniprot_id][1]
            if uniprot_id!=old_id:
                i=1
            else:
                i+=1
            poly_name = uniprot_id+"_p_"+str(i)
            poly_num = i
            poly_start = int(rowsplit[1])
            poly_rel_start = poly_start/seq_len
            poly_end = int(rowsplit[2])
            poly_rel_end = poly_end/seq_len
            poly_parts = rowsplit[3]
            poly_part = rowsplit[3].split("/")
            poly_part1 = poly_part[0]
            poly_part1_tp = resources.get_aa_type(poly_part1)
            if (tp==2):
                poly_part2 = poly_part[1]
                poly_part2_tp = resources.get_aa_type(poly_part2)
                if (poly_part1_tp==poly_part2_tp):
                    poly_pair_tp = poly_part1_tp
                else:
                    poly_pair_tp = "Both"
                poly_aa = rowsplit[4]
            else:
                poly_part2=""
                poly_part2_tp=""
                poly_pair_tp=poly_part1_tp
                poly_aa = rowsplit[6]
            poly_size = len(poly_aa)
            idr_group_tots, idr_group_props = resources.get_AAcomposition(poly_aa, 2)
            poly_lst.append([uniprot_id, poly_name, poly_num, poly_start, poly_end, 
                             poly_rel_start, poly_rel_end, poly_parts, poly_part1, 
                             poly_part2, poly_part1_tp, poly_part2_tp, poly_pair_tp, poly_aa, 
                             poly_size, seq_len, seq_desc, seq_aa]+idr_group_tots+idr_group_props)
            old_id = uniprot_id
    return poly_lst


### POLY PAIR analysis ###
def get_poly_mask(col1, col2, col3):
    ''' Generates a mask to facilitate the grouping of categories of  PolyXYs. '''
    return col3.replace(col1, 'X').replace(col2, 'Y')


def get_poly_patt(col1):
    ''' Counts the successive residues to group palindromic or different kinds of pairs. '''
    c_ant = ""
    patt = ""
    for c in col1:
        if c_ant!=c:
            if c_ant!="":
                if patt=="":
                    patt=patt+str(ct)
                else:
                    patt=patt+"."+str(ct)
            ct=1
        else:
            ct+=1
        c_ant=c
    patt=patt+"."+str(ct)
    return patt


def get_poly_cat(col1):
    ''' Create different types of categories depending on how the repeat is organized. '''
    s_counts = [int(x) for x in col1.split(".")]
    s_twos = sum(True if x>=2 else False for x in s_counts)
    s_ones = all(x == 1 for x in s_counts)
    s_ones_plus = all(x == 1 for x in s_counts[:-1]) or all(x == 1 for x in s_counts[1:])
    s_pall = False
    if (len(s_counts)%2!=0):
        midd = math.ceil(len(s_counts)/2)
        rev = s_counts[midd:]
        rev = rev[::-1]
        s_pall = s_counts[:midd-1]==rev
    
    if s_ones:
        poly_cat = "(XY)n"
    elif s_ones_plus:
        poly_cat = "(XY)n+"
    elif len(s_counts)==2:
        poly_cat = "[X]n[Y]z"
    elif s_pall:
        poly_cat = "n[XY]n"
    else:
        poly_cat = "[XY]n"
    return poly_cat


# df_poly_details = df_poly_details.drop(['poly_pairs', 'poly_pairs1', 'poly_pairs2', 'poly_pairs1_int', 'poly_pairs2_int', 'poly_mask', 'poly_patt', 'poly_cat'], axis=1)
def generate_poly_groups(df_poly_details):
    ''' Prepare and add the categories of poly to dataframe. '''
    df_poly_details['poly_pairs'] = df_poly_details.apply(lambda x: resources.get_pairs(x.poly_part1, x.poly_part2), axis=1)
    df_poly_details[['poly_pairs1', 'poly_pairs2']] = df_poly_details['poly_pairs'].str.split('|', 1, expand=True)
    df_poly_details['poly_pairs1_int'] = df_poly_details['poly_pairs1'].apply(lambda x: resources.AA_CODE_DICT[x])
    df_poly_details['poly_pairs2_int'] = df_poly_details['poly_pairs2'].apply(lambda x: resources.AA_CODE_DICT[x])
    df_poly_details['poly_mask'] = df_poly_details.apply(lambda x: get_poly_mask(x.poly_pairs1, x.poly_pairs2, x.poly_aa), axis=1)
    df_poly_details['poly_patt'] = df_poly_details.apply(lambda x: get_poly_patt(x.poly_mask), axis=1)
    df_poly_details['poly_cat'] = df_poly_details.apply(lambda x: get_poly_cat(x.poly_patt), axis=1)
    return df_poly_details

### IDR COVERAGE by POLY ###
def extract_poly_idr(df_idr_details, df_poly, cutoff=.6, min_sz=6):

    idrs_pos = df_idr_details.loc[:, ['seq_name', 'idr_name', 'idr_size', 'idr_start', 'idr_end']]
    idrs_pos = idrs_pos.rename(columns={"over_sz": "over_idr_sz"})
    merged_pos = pd.merge(df_poly, idrs_pos, how="left", on="seq_name")
    merged_pos = merged_pos.sort_values(by=['idr_name', 'poly_name'])
    # Here instead of using the IDR real start and ends I use the ones calculated
    # based on the alignment with the PDB structure (avoid wrong overlaps)
    starts_ends = merged_pos.loc[:, ['idr_start', 'idr_end', 
                                     'poly_start', 'poly_end']].values
    poly_sizes = merged_pos['poly_size'].values
   
    pos_over, sz_over = cross.define_overlaps(starts_ends)
    orig_over_perc = cross.get_overlap_size(sz_over, poly_sizes)
    _, long_bool = cross.extract_short_overlaps(sz_over, orig_over_perc, cutoff, min_sz)
    merged_sel = merged_pos.copy()
    merged_sel = merged_sel.loc[long_bool, :]
    
    # Adding the just calculated over_sz for the poly region to use to determine
    # the SS region (gaps will be accounted there)
    sz_over = sz_over[long_bool]
    orig_over_perc = orig_over_perc[long_bool]
    pos_over = pos_over[long_bool, :]
    merged_sel.loc[:, 'over_sz_idr'] = sz_over
    merged_sel.loc[:, 'over_perc_idr'] = orig_over_perc
    merged_sel.loc[:, 'poly_start_idr'] = (pos_over[:,0]).astype(int)
    merged_sel.loc[:, 'poly_end_idr'] = (pos_over[:,1]).astype(int)
    merged_sel.loc[:, 'region_idr'] = merged_sel.apply(lambda x: x['seq_aa'][int(x['poly_start_idr']-1):int(x['poly_end_idr'])], axis=1)
    merged_sel[['tot_polar_idr','tot_non_polar_idr']] = merged_sel.apply(lambda x: resources.get_AAcomposition(x['region_idr']), axis=1)
    merged_sel = merged_sel.drop(['region_idr'], axis=1)
    
    # Adding the polyXYs without PDB overlaps. I use the not in from the list
    # of selected because the merge caused duplications
    IDs_seq = pd.unique(merged_pos.loc[long_bool, 'poly_name']).tolist()
    df_poly_notSel = df_poly.loc[~df_poly['poly_name'].isin(IDs_seq), :]
    df_poly_new = pd.concat([merged_sel, df_poly_notSel])
    
    return df_poly_new


def count_seq_pairs(sequence, k_mer=2, acc_pairs=True):
    """Generate the slide for each integer for each amino-acid"""
    seq_int_slide = []
    sequence = [aa for aa in sequence]
    seq = np.array(sequence)
    n = len(sequence)
    xcol = np.arange(n-k_mer+1)
    inds_col = np.tile(np.arange(k_mer), (n-k_mer+1,1))
    inds_row = np.tile(xcol, (k_mer,1)).transpose()
    seq_int_slide_ori = seq[inds_col + inds_row]
    #    inds2 = xcol.reshape(int((n-1)/2), 2)
    #    seq_int_slide_ori = seq[inds2]
    seq_int_slide = seq_int_slide_ori.copy()
    # I'll accumulate the inverse pairs (XY and YX) for now, but make it flexible
    # if we decide not to accumulate in te future.
    if acc_pairs:
        check_inv = seq_int_slide_ori[:,0]>seq_int_slide_ori[:,1]
        seq_int_slide[check_inv, 0] = seq_int_slide_ori[check_inv,1]
        seq_int_slide[check_inv, 1] = seq_int_slide_ori[check_inv,0]        
    #check_dif = seq_int_slide[:,0]!=seq_int_slide[:,1]
    #seq_int_slide = seq_int_slide[check_dif, :]
    pairs_rep = np.apply_along_axis(''.join, 1, seq_int_slide)
    pairs = pd.DataFrame(np.transpose(np.array(np.unique(pairs_rep, return_counts=True))),
                         columns=['poly_pair_simp', 'cnt'])
    pairs = pairs.set_index('poly_pair_simp')
    pairs = pairs.astype(int)
    return pairs


def count_pairs(polyXYfas_path, k_mer=2, acc_pairs=True):
    seq_dict = cross.load_seqs(polyXYfas_path)
    i=0
    for vals in seq_dict.values():
        i+=1
        pairs = count_seq_pairs(vals[0], k_mer, acc_pairs)
        if i==1:
            df_pairs = pairs.copy()
        else:
            df_pairs = df_pairs.merge(pairs, how="outer", on="poly_pair_simp").fillna(0)
            #if "QX" in df_pairs.index:
            #    break
            df_pairs = pd.DataFrame(df_pairs.sum(axis=1), columns=['cnt'])
        df_pairs = df_pairs.sort_index()
    df_pairs['prop'] = df_pairs['cnt']/int(df_pairs.sum(axis=0))
    #df_pairs = df_pairs.sort_index()
    df_pairs = df_pairs.sort_values(by='prop')
    df_pairs = df_pairs.reset_index()
    df_pairs['poly_pairs1_int'] = df_pairs['poly_pair_simp'].apply(lambda x: AA_CODE_DICT[x[0]])
    df_pairs['poly_pairs2_int'] = df_pairs['poly_pair_simp'].apply(lambda x: AA_CODE_DICT[x[1]])
    return df_pairs.iloc[:, [0,-2,-1,1,2]]

def count_aas(seqs_path):
    seq_dict = cross.load_seqs(seqs_path)
    i=0
    cnt_mat = np.empty([len(seq_dict), len(AA_CODE_LIST)], dtype=int)
    for k, vals in seq_dict.items():
        #print(k + " - " + str(i))
        chars = np.array([AA_CODE_LIST.index(r) for r in vals[0]])
        cnt_mat[i] = np.bincount(chars, minlength=len(AA_CODE_LIST))
        i += 1
    tot = np.sum(cnt_mat, 0)
    tot_frac = tot/np.sum(tot)
    tot = np.stack([AA_CODE_LIST, tot, tot_frac], 1)
    tot = np.delete(tot, 0, 0)
    df_aas = pd.DataFrame(tot, columns=['poly_pairs','cnt', 'prop'])
    df_aas = df_aas.astype({'cnt': 'int32', 'prop': 'float32'})
    return (df_aas)
        

### POLY - IDR - PDB relations ###
def extract_poly_idr_pdb(df_idr_details, df_poly, cutoff=.6, min_sz=4):

    # Filtering the PDB important columns from the IDR perspective
    df_idr_details = df_idr_details.loc[df_idr_details['idr_type']=='predicted', :]
    idrs_pos = df_idr_details.loc[:, ['idr_name', 'ss_final_over_sz', 'pdb_name', 'pdb_evalue',
       'bitscore', 'internal_id', 'internal_id2', 'seq_start', 'seq_end', 'pdb_start', 
       'pdb_end', 'hit_id', 'hsp_id', 'removed_on_s', 'removed_on', 'kept_on', 
       'idr_seq_rel_start', 'idr_seq_rel_end', 'idr_final_seq_align_start', 
       'idr_final_seq_align_end']]
    idrs_pos = idrs_pos.rename(columns={"ss_final_over_sz": "ss_over_idr_sz"})
    merged_pos = pd.merge(df_poly, idrs_pos, how="left", on="idr_name")
    merged_pos = merged_pos.sort_values(by=['idr_name', 'poly_name'])
    # Here instead of using the IDR real start and ends I use the ones calculated
    # based on the alignment with the PDB structure (avoid wrong overlaps)
    starts_ends = merged_pos.loc[:, ['idr_seq_rel_start', 'idr_seq_rel_end', 
                                     'poly_start', 'poly_end']].values
    poly_sizes = merged_pos['poly_size'].values
   
    pos_over, sz_over = cross.define_overlaps(starts_ends)
    orig_over_perc = cross.get_overlap_size(sz_over, poly_sizes)
    _, long_bool = cross.extract_short_overlaps(sz_over, orig_over_perc, cutoff, min_sz)
    merged_sel = merged_pos.copy()
    merged_sel = merged_sel.loc[long_bool, :]
    
    # Adding the just calculated over_sz for the poly region to use to determine
    # the SS region (gaps will be accounted there)
    sz_over = sz_over[long_bool]
    orig_over_perc = orig_over_perc[long_bool]
    merged_sel.loc[:, 'over_sz'] = sz_over
    merged_sel.loc[:, 'over_perc'] = orig_over_perc
    
    # Adding the polyXYs without PDB overlaps. I use the not in from the list
    # of selected because the merge caused duplications
    IDs_seq = pd.unique(merged_pos.loc[long_bool, 'poly_name']).tolist()  
    df_poly_notSel = df_poly.loc[~df_poly['poly_name'].isin(IDs_seq), :]
    df_poly_new = pd.concat([merged_sel, df_poly_notSel])
    
    return df_poly_new


def merge_all_coverages(sel_ids, sz_over):
    _, un_idx, rep_idx = np.unique(sel_ids, return_index=True, return_inverse=True)
    new_sz_over = np.bincount(un_idx[rep_idx], weights=sz_over)
    new_sz_over = new_sz_over[un_idx]
    return new_sz_over, un_idx, rep_idx


def merge_poly_coverage(df_poly_details, cols, df_idr_details, comp="idr", par="idr_name"):
    all_cols = ['poly_name', 'poly_pairs', 'poly_pair_tp', par]+cols
    df_poly_sel = (df_poly_details.loc[~df_poly_details[par].isna(), all_cols].
                   sort_values(by=['idr_name', cols[-2], cols[-1]]))
    idr_ids = df_poly_sel.loc[:,'idr_name'].values
    un_idrs, un_ids_idrs, inv_ids, counts_poly = np.unique(idr_ids, return_counts=True, 
                                                           return_index=True, 
                                                           return_inverse=True)
    poly_ids = df_poly_sel.loc[:,'poly_name'].values
    poly_pairs = df_poly_sel.loc[:,'poly_pairs'].values
    poly_polar_cat = df_poly_sel.loc[:,'poly_pair_tp'].values
    poly_polar = df_poly_sel.loc[:,cols[-4]].values
    poly_non_polar = df_poly_sel.loc[:,cols[-3]].values
    poly_coords = ["|".join(i) for i in df_poly_sel.loc[:,cols[-2:]].values.astype(str)]
    if (comp=="pdb"):
      poly_coords_ss = ["|".join(i) for i in df_poly_sel.loc[:,cols[-6:-4]].values.astype(str)]
      poly_coords_ss_d = dict()
    un_ids_pairs = np.append(un_ids_idrs[1:], len(poly_ids))
    polyID_d, poly_pairs_d, poly_polar_cat_d, poly_coords_d = dict(), dict(), dict(), dict()
    poly_polar_d, poly_non_polar_d = dict(), dict()
    for i,j,k in zip(un_idrs, un_ids_idrs, un_ids_pairs):
        polyID_d[i] = list(poly_ids[j:k])
        poly_pairs_d[i] = list(poly_pairs[j:k])
        poly_polar_cat_d[i] = list(poly_polar_cat[j:k])
        poly_coords_d[i] = list(poly_coords[j:k])
        poly_polar_d[i] = sum(poly_polar[j:k])
        poly_non_polar_d[i] = sum(poly_non_polar[j:k])
        if (comp=="pdb"):
            poly_coords_ss_d[i] = list(poly_coords_ss[j:k])
        
    sz_over, un_idx, rep_idx = merge_all_coverages(idr_ids, df_poly_sel[cols[0]].values)
    idr_sizes = df_poly_sel.loc[:,cols[1]].values
    idr_sizes = idr_sizes[un_idx]
    over_idr = cross.get_overlap_size(sz_over, idr_sizes)
    
    df_over_idr = pd.DataFrame(np.hstack((un_idrs[:, None], over_idr[:, None],
                                          sz_over[:, None], counts_poly[:, None])), 
                               columns=["idr_name", "poly_over_"+comp, "poly_sum_sz_"+comp, "poly_cnt_"+comp])
    df_over_idr["poly_over_"+comp] = df_over_idr["poly_over_"+comp].astype(float)
    #df_over_idr["poly_sum_sz_"+comp] = df_over_idr["poly_sum_sz_"+comp].astype(int)
    df_over_idr["polyXYs_"+comp] = df_over_idr["idr_name"].apply(lambda x: ", ".join(polyID_d[x]))
    df_over_idr["poly_pairs_"+comp] = df_over_idr["idr_name"].apply(lambda x: ", ".join(poly_pairs_d[x]))
    df_over_idr["poly_polarity_"+comp] = df_over_idr["idr_name"].apply(lambda x: ", ".join(poly_polar_cat_d[x]))
    df_over_idr["poly_tot_polar_"+comp] = df_over_idr["idr_name"].apply(lambda x: poly_polar_d[x])
    df_over_idr["poly_tot_non_polar_"+comp] = df_over_idr["idr_name"].apply(lambda x: poly_non_polar_d[x])
    df_over_idr["poly_coords_"+comp] = df_over_idr["idr_name"].apply(lambda x: ", ".join(poly_coords_d[x]))
    if (comp=="pdb"):
        df_over_idr["poly_coords_ss"] = df_over_idr["idr_name"].apply(lambda x: ", ".join(poly_coords_ss_d[x]))
    df_idr_details = pd.merge(df_idr_details, df_over_idr, how="left", on="idr_name")
    return df_idr_details


def merge_poly_ss_coverage(df_2D_details, df_poly_details, df_idr_details):
    df_ss_sel = df_2D_details.loc[:, ['poly_name', 'poly_ss']]
    df_poly_sel = df_poly_details.loc[:, ['poly_name', 'idr_name', 'poly_start']]
    df_ss_sel = df_ss_sel.merge(df_poly_sel, how='left', on='poly_name')
    df_ss_sel = df_ss_sel.sort_values(by=['idr_name', 'poly_start'])
    
    poly_regions = df_ss_sel.loc[:,'poly_ss'].values
    poly_regions_d = dict()
    un_idrs, un_ids_idrs, inv_ids, counts_poly = np.unique(df_ss_sel.idr_name.values, return_counts=True, 
                                                           return_index=True, return_inverse=True)
    un_ids_pairs = np.append(un_ids_idrs[1:], len(poly_regions))
    for i,j,k in zip(un_idrs, un_ids_idrs, un_ids_pairs):
        poly_regions_d[i] = list(poly_regions[j:k])
    
    df_ss = pd.DataFrame(df_ss_sel.loc[df_ss_sel['idr_name'].isin(un_idrs), 'idr_name'].drop_duplicates())
    df_ss["poly_ss_region_un"] = df_ss['idr_name'].apply(lambda x: "".join(poly_regions_d[x]))
    df_2D = cross.prepare_ss_counts(df_ss, "", "pdb")
    df_idr_details = df_idr_details.merge(df_2D, how="left", on="idr_name")
    
    return df_idr_details


def generate_poly_mask(idr_name, seq_aa, poly_coords_idr, idr_start, idr_end, side):
    #print(idr_name)
    sub=0
    if (side!="_ss"):
        sub=1
    poly_coords = poly_coords_idr.split(', ')
    for i in range(len(poly_coords)):
        pair = poly_coords[i].split('|')
        start = int(float(pair[0]))
        end = int(float(pair[1]))
        seq_aa = seq_aa[:start-sub]+'|'*(end-start+sub)+seq_aa[end:]
    seq_aa = '|'*(int(idr_start)-sub)+seq_aa[int(idr_start)-sub:int(idr_end)]+'|'*(len(seq_aa)-int(idr_end))
    return seq_aa


def extract_region_noPoly(df_idr_details, cols, side=""):
    #cols = ['seq_aa', 'poly_coords_idr', 'idr_start', 'idr_end']
    if (side!=""):
        side = "_"+side
    df_sel = df_idr_details.loc[~df_idr_details[cols[1]].isna(), ['idr_name']+cols]
    df_sel['idr_poly_mask'+side] = (df_sel.apply(lambda x: 
                                            generate_poly_mask(x['idr_name'],
                                                               x[cols[0]],
                                                               x[cols[1]],
                                                               x[cols[2]],
                                                               x[cols[3]], side), axis=1))
    df_sel['idr_noPoly'+side] = df_sel['idr_poly_mask'+side].str.replace('|', '')
    if (side=="_ss"):
        df_ss = df_sel.loc[:, ['idr_name', 'idr_noPoly'+side]]
        df_2D = cross.prepare_ss_counts(df_ss, "", "noPDB")
        df_idr_details = df_idr_details.merge(df_2D, how="left", on="idr_name")
    else:
        df_sel[['tot_polar_noPoly'+side,'tot_non_polar_noPoly'+side]] = df_sel.apply(lambda x: resources.get_AAcomposition(x['idr_noPoly'+side]), axis=1)
    df_sel = df_sel.drop(cols, axis=1)
    df_idr_details = pd.merge(df_idr_details, df_sel, how="left", on="idr_name")
    return df_idr_details


def get_align_aa_counts(align_reg, pair1, pair2):
    ''' Counts how many occurrences of each AA appear after the alignment with PDB.'''
    c = collections.Counter(align_reg)
    return pd.Series((c[pair1], c[pair2]))


### MAIN POLY SESSION
def run_poly(seqs_path, polyXY_path, n_aa, cutoff, min_size):
    fasta_data,_,_ = resources.read_fasta(seqs_path)
    poly_lst = extract_polyXY_features(polyXY_path, fasta_data, n_aa)
    colnames = ["seq_name", "poly_name", "poly_num", "poly_start", "poly_end", 
                "poly_rel_start", "poly_rel_end", "poly_parts", "poly_part1", 
                "poly_part2", "poly_part1_tp", "poly_part2_tp", "poly_pair_tp", 
                "poly_aa", "poly_size", "seq_len", "seq_desc", "seq_aa", 
                "tot_polar", "tot_non_polar", "prop_polar", "prop_non_polar"]
    pd_poly = pd.DataFrame(poly_lst, columns=colnames)
    pd_poly = pd_poly.sort_values(by=['seq_name', 'poly_start', 'poly_end'])
    
    if (n_aa==1):
        #Use it when processing just PolyX. This data will not be used
        source = "polyx"
        pd_poly = pd_poly.drop(columns={"prop_polar", "prop_non_polar", "poly_part1", "poly_part2", "poly_part1_tp", "poly_part2_tp"})
    else:
        # Don't run this for polyX. Makes sense when we have more than 2 different AAs
        source = "polyxy"
        pd_poly = generate_poly_groups(pd_poly)
    
    # Annotating the polyX/Ys inside IDRs
    #pd_poly = extract_poly_idr(df_idr_details, pd_poly, cutoff, min_size)
    polycsv_path = resources.gen_filename(polyXY_path, source, "details", "csv")
    pd_poly.to_csv(polycsv_path, index=False)
    
    return (polycsv_path)


### MAIN IDR POLY SESSION
def run_idr_poly(df_idr_details, df_poly, polyXYfas_path, pdb_mask_path, dssp_path, blast_over_path, path_sel_coords, ssSeq_path, polyidr_path, polyss_path):
    #df_idr_details = pd.read_csv(idrcsv_path)
    # Poly over IDR regardless of PDB structure
    
    df_poly_details[['cnt_pair1','cnt_pair2']] = df_poly_details.apply(lambda x: get_align_aa_counts(x['poly_aa'], x['poly_pairs1'], x['poly_pairs2']), axis=1)
        
    # Getting IDR side of the overlaps (list of PolyXYs and coverage with PDB)
    df_poly_details = extract_poly_idr_pdb(df_idr_details, df_poly_details)
        
    # get the details of the PDB and SEQ alignment
    df_poly_details = cross.append_pdb_details(pdb_mask_path, df_poly_details)
    # Now we can calculate the PolyXY ss_region using the same function
    df_poly_details, df_2D_details, _ = cross.get_dssp_idr(dssp_path, blast_over_path, df_poly_details, 50)
    df_poly_details = cross.correct_pdb_coords(path_sel_coords, df_poly_details, ssSeq_path, polyidr_path)
    
    df_poly_details_simp = df_poly_details.loc[~df_poly_details["pdb_name"].isna(), ["poly_name", "poly_pairs1", "poly_pairs2"]]
    df_2D_details = pd.merge(df_2D_details, df_poly_details_simp, how="left", on="poly_name")
    df_2D_details[['cnt_pair1_align','cnt_pair2_align']] = df_2D_details.apply(lambda x: get_align_aa_counts(x['poly_seq'], x['poly_pairs1'], x['poly_pairs2']), axis=1)
    df_2D_details = df_2D_details.drop({"poly_pairs1", "poly_pairs2"}, axis=1)
    df_2D_details.to_csv(polyss_path, index=False)
    return df_poly_details, df_2D_details


def run_idr_poly_global(df_idr_details, df_poly):
    # Getting IDR coverage area based on all PolyXYs over the IDR
    cols = ['over_sz_idr', 'idr_size', 'tot_polar_idr', 'tot_non_polar_idr', 
            'poly_start_idr', 'poly_end_idr']
    df_idr_details = merge_poly_coverage(df_poly_details, cols, df_idr_details)
    
    cols = ['seq_aa', 'poly_coords_idr', 'idr_start', 'idr_end']
    df_idr_details = extract_region_noPoly(df_idr_details, cols)
    
    # Getting IDR side of the overlaps (list of PolyXYs and coverage with PDB)
    comp = "pdb"
    cols = ['ss_over_sz', 'ss_over_idr_sz', 'idr_name', 
            'poly_seq_align_start', 'poly_seq_align_end','tot_polar_pdb', 
            'tot_non_polar_pdb', 'poly_seq_rel_start', 'poly_seq_rel_end']
    df_idr_details = merge_poly_coverage(df_poly_details, cols, df_idr_details, comp, 'pdb_name')
    cols = ['seq_aa', 'poly_coords_pdb', 'idr_seq_rel_start', 'idr_seq_rel_end']
    df_idr_details = merge_poly_ss_coverage(df_2D_details, df_poly_details, df_idr_details)
    df_idr_details = extract_region_noPoly(df_idr_details, cols, 'pdb')
    
    # Extracting the SS region for the IDR region not covered by PolyXYs
    cols = ['final_align_ss', 'poly_coords_ss', 'idr_seq_align_start', 'idr_seq_align_end']
    df_idr_details = extract_region_noPoly(df_idr_details, cols, 'ss')
    
    cols = list(df_idr_details.columns[[i for i,x in enumerate(df_idr_details.columns) if x == "poly_over_idr"][0]:])
    df_idr_details = df_idr_details.loc[~df_idr_details['poly_over_idr'].isna(), ["idr_name"]+cols]
    #df_idr_details.to_csv(idrpoly_path, index=False)
    
    # Counting the Poly pairs (this doesn't change)
    df_pairs_poly = count_pairs(seqs_path, 2, True)
    # Pablo suggested me to count all to check for possible bias in sequences with repeated regions
    # As expected the difference is thin. I'll check with both.
    df_pairs_all = count_pairs(seqs_path, 2, True)
    df_pairs_all.to_csv('analysis/data_pairs_all.csv', index=False)
    ssprop_lst, ss_counts, ss_props = cross.count_dssp_seqs(dssp_path, df_poly_details)
    
    df_aa_all = count_aas(seqs_path)
    df_aa_all.to_csv('analysis/data_aa_all.csv', index=False)
    

def run_properties():
    _ = get_idr_properties.get_cider_props(df_poly_details, polypar_path, pref)

pref='poly'
#df_poly = pd.read_csv(polycsv_path)
#df_poly_new = pd.read_csv('analysis_masked/data_noSeqs_'+pref+'.csv')
#df_poly_details = pd.read_csv(polycsv_path)
#df_poly_details = pd.read_csv(polyidr_path)

#file_names=['evalues_'+pref+'.pickle', 'idr_sizes_'+pref+'.pickle', 
#                'unique_ids_'+pref+'.pickle', 'gen_info_'+pref+'.pickle', 
#                'blast_over_'+pref+'.pickle']
