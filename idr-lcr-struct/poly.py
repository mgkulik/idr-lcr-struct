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
import idrPdb as cross
        

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
            # PolyX2 has a complete different order of columns now
            rowsplit = line[1].rstrip("\n").split("\t")
            if tp==1:
                rowsplit = [rowsplit[5], rowsplit[0], rowsplit[1], rowsplit[2], rowsplit[3], rowsplit[4], rowsplit[6]]
            uniprot_id = rowsplit[0].split("|")[1]
            seq_val = fasta_data[uniprot_id]
            seq_aa = str(seq_val.seq)
            seq_len = len(seq_aa)
            seq_desc = seq_val.description
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
            if (poly_part1!="X"):
                poly_part1_tp = resources.get_aa_type(poly_part1)
            if (tp==2):
                poly_part2 = poly_part[1]
                if (poly_part2!="X"):
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
            if (poly_part1!="X") and (poly_part2!="X"):
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
    df_poly_details[['cnt_pair1','cnt_pair2']] = df_poly_details.apply(lambda x: get_align_aa_counts(x['poly_aa'], x['poly_pairs1'], x['poly_pairs2']), axis=1)
    return df_poly_details


### IDR COVERAGE by POLY ###
def extract_poly_idr(idrcsv_path, df_poly, cutoff=.6, min_size=6):
    ''' Add the information of all IDRs to the poly annotations. '''

    df_idr_details = pd.read_csv(idrcsv_path)
    idrs_pos = df_idr_details.loc[:, ['seq_name', 'idr_name', 'idr_size', 'idr_start', 'idr_end']]
    merged_pos = pd.merge(df_poly, idrs_pos, how="left", on="seq_name")
    merged_pos = merged_pos.sort_values(by=['idr_name', 'poly_name'])
    # Here instead of using the IDR real start and ends I use the ones calculated
    # based on the alignment with the PDB structure (avoid wrong overlaps)
    starts_ends = merged_pos.loc[:, ['idr_start', 'idr_end', 
                                     'poly_start', 'poly_end']].values
    poly_sizes = merged_pos['poly_size'].values
   
    pos_over, sz_over = cross.define_overlaps(starts_ends)
    orig_over_perc = cross.get_overlap_size(sz_over, poly_sizes)
    _, long_bool = cross.extract_short_overlaps(sz_over, orig_over_perc, cutoff, min_size)
    merged_sel = merged_pos.copy()
    merged_sel = merged_sel.loc[long_bool, :]
    
    # Adding the just calculated over_sz for the poly region to use to determine
    # the SS region (gaps will be counted later)
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


### POLY - IDR - PDB relations ###
def extract_poly_idr_pdb(df_idr_details, df_poly, cutoff=.6, min_size=4):
    ''' Filter all cases that overlap IDRs and PDB sequences. Important fields
    from the IDRs are collected. '''
    
    # Filtering the PDB important columns from the IDR perspective
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
    _, long_bool = cross.extract_short_overlaps(sz_over, orig_over_perc, cutoff, min_size)
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
    df_poly_new["pdb_id"] = df_poly_new.pdb_name.str.split("_", n=1, expand=True)[0]
    
    return df_poly_new


def get_align_aa_counts(align_reg, pair1, pair2):
    ''' Counts how many occurrences of each AA appear after the alignment with PDB.'''
    c = collections.Counter(align_reg)
    return pd.Series((c[pair1], c[pair2]))

def final_2Ddata(df_poly_details, df_2D_details):
    df_poly_details_simp = df_poly_details.loc[~df_poly_details["pdb_name"].isna(), ["poly_name", "poly_pairs1", "poly_pairs2"]]
    df_2D_details = pd.merge(df_2D_details, df_poly_details_simp, how="left", on="poly_name")
    df_2D_details[['cnt_pair1_align','cnt_pair2_align']] = df_2D_details.apply(lambda x: get_align_aa_counts(x['poly_seq'], x['poly_pairs1'], x['poly_pairs2']), axis=1)
    df_2D_details = df_2D_details.drop({"poly_pairs1", "poly_pairs2"}, axis=1)
    return df_2D_details


def subset_fasta(fasta_data, df_poly_details):
    ''' Creating a fasta with the subset of uniprot ids of interest. '''
    seq_lst = []
    un_seq_names = df_poly_details.loc[~df_poly_details["pdb_name"].isna(), "seq_name"]
    inter_seqs = np.intersect1d(un_seq_names, list(fasta_data.keys())).tolist()
    for i in range(0, len(inter_seqs)):
        seq_lst.append(fasta_data[inter_seqs[i]])
    return (seq_lst)


### MAIN POLY SESSION
def main_poly(seqs_path, polyXY_path, idrcsv_path, source, n_aa, cutoff, min_size):
    ''' Collecting poly info and their relation to IDRs. '''
    fasta_data,_,_ = resources.read_fasta(seqs_path)
    poly_lst = extract_polyXY_features(polyXY_path, fasta_data, n_aa)
    colnames = ["seq_name", "poly_name", "poly_num", "poly_start", "poly_end", 
                "poly_rel_start", "poly_rel_end", "poly_parts", "poly_part1", 
                "poly_part2", "poly_part1_tp", "poly_part2_tp", "poly_pair_tp", 
                "poly_aa", "poly_size", "seq_len", "seq_desc", "seq_aa", 
                "tot_polar", "tot_non_polar", "prop_polar", "prop_non_polar"]
    pd_poly = pd.DataFrame(poly_lst, columns=colnames)
    pd_poly['proteome_poly_id'] = resources.extract_proteome(seqs_path)
    pd_poly = pd_poly.sort_values(by=['seq_name', 'poly_start', 'poly_end'])
    
    if (source == "polyx"):
        #Use it when processing just PolyX. This data will not be used
        pd_poly = pd_poly.drop(columns={"prop_polar", "prop_non_polar", "poly_part1", "poly_part2", "poly_part1_tp", "poly_part2_tp"})
        
    else:
        # Don't run this for polyX. Makes sense when we have more than 2 different AAs
        pd_poly = generate_poly_groups(pd_poly)
    
    # Annotating the polyX/Ys inside IDRs
    pd_poly = extract_poly_idr(idrcsv_path, pd_poly, cutoff, min_size)
    polycsv_path = resources.gen_filename(polyXY_path, source, "details", "csv")
    pd_poly.to_csv(polycsv_path, index=False)
    
    return (polycsv_path)


### MAIN IDR POLY SESSION
def main_poly_pdb(idr_all_path, poly_details_path, pdb_mask_path, dssp_path, file_names, path_coords, path_fasta, source, cutoff, min_size):
    ''' Now crossing Poly with PDBs. '''
    
    basis_path = resources.get_dir(idr_all_path)+"/"+resources.get_filename(idr_all_path)+"_"
    
    if len(file_names)==5:
        file_names.insert(0, 'starts_ends')

    df_idr_details = pd.read_csv(idr_all_path, low_memory=False)
    df_poly_details = pd.read_csv(poly_details_path, low_memory=False)
    # Getting IDR side of the overlaps (list of PolyXYs and coverage with PDB)
    df_poly_details = extract_poly_idr_pdb(df_idr_details, df_poly_details, cutoff, min_size)
        
    # get the details of the PDB and SEQ alignment
    df_poly_details = cross.append_pdb_details(pdb_mask_path, df_poly_details)
    # Now we can calculate the PolyXY ss_region using the same functions. The filtering step is not required
    df_poly_details, df_2D_details, _ = cross.get_dssp_idr(dssp_path, basis_path+file_names[-1]+".pickle", df_poly_details, basis_path, "poly", 50, source)
    poly_all_path = "data_all_"+source+".csv"
    _ = cross.mark_no_homologs(basis_path, df_poly_details, poly_all_path)
    df_pdb_coords = pd.read_csv(path_coords)
    cross.correct_pdb_coords(df_pdb_coords, df_poly_details, poly_all_path, "poly")
    
    if (source=="polyxy"):
        df_2D_details = final_2Ddata(df_poly_details, df_2D_details)
    polyss_path = basis_path+"data_ss_"+source+".csv"
    df_2D_details.to_csv(polyss_path, index=False)
    
    # Generate a subset of the fasta with only the sequences overlaping the 3 sets
    fasta_data,_,_ = resources.read_fasta(path_fasta)
    seq_lst = subset_fasta(fasta_data, df_poly_details)
    fastaout = basis_path+"idr_"+source+"_pdb.fasta"
    resources.save_fastas(seq_lst, fastaout)
    
    return basis_path+poly_all_path, polyss_path