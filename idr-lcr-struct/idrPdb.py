#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 14 16:38:58 2020

@author: magoncal
"""

from Bio import SeqIO, SearchIO
from Bio.SubsMat.MatrixInfo import blosum62 as blosum
from collections import defaultdict
import pandas as pd
import numpy as np
import lxml.etree as ET

import pickle
import time
import os
import re
import math

import resources

pd.options.display.max_columns = 50

STRUC2D_DICT = {'H':0,'G':0,'I':0,'B':1,'E':1,'T':2,'S':2,'C':2,'X':3,' ':4,'-':5, '|':6}
STRUC2D_GRP = ['alpha-helix', 'beta-sheet','coil','unmodeled','unfolded', 'gaps', 'not aligned']

AA_CODE_LIST = ['?','A','R','N','D','C','Q','E','G','H','I','L','K','M','F','P','S','T','W','Y','V','U','O','B','Z','X']
AA_CODE_DICT = {'A':1,'C':2,'D':3,'E':4,'F':5,'G':6,'H':7,'I':8,'K':9,'L':10,'M':11,'N':12,'P':13,'Q':14,'R':15,'S':16,'T':17,'V':18,'W':19,'Y': 20,'U':21,'O':22,'B':23,'Z':24,'X':25}
AA_GROUPS = [[19,8,16,17,5,3,6,22,12,2,9,4,7], [14,18,1,15,20,11,10,13,21]]
AA_GROUPS_NAMES = ['Polar', 'Non-Polar']

blosum.update(((b,a),val) for (a,b),val in list(blosum.items()))
# Added U pairs as zero. As it is uncommon this should not be a problem. 
blosum["U", "T"] = 0
blosum["U", "C"] = 0

##### PART 1, EXTRACTING THE DATA FROM BASE FILES AND SAVING TO DISK #####

def extract_blast_pairs(xml_path, pickle_sz, sep1, sep2):
    ''' Gets the blastP xml results and generates pickle file with the 
    overlaps data. '''
    
    pdb_filename = resources.gen_filename(xml_path, "blast", "", "pickle")
    
    hsp_cols = ['Hsp_num', 'Hsp_bit-score', 'Hsp_evalue', 
            'Hsp_query-from', 'Hsp_query-to', 'Hsp_hit-from', 'Hsp_hit-to',
            'Hsp_identity', 'Hsp_gaps', 'Hsp_qseq', 'Hsp_hseq', 'Hsp_midline']
    
    i=0
    mode =""
    dt_pdb = dict()
    start_time = time.time()
    for event, elem in ET.iterparse(xml_path, events=('end',)):
        if elem.tag == "Iteration_query-def":
            query_vals = list()
            query_vals.append(elem.text.split(sep1)[1])
            elem.clear()
        elif elem.tag == "Hit":
            hits_vals = list()
            for hit_child in elem:
                if hit_child.tag == "Hit_num":
                    hits_vals.append(hit_child.text)
                elif hit_child.tag == "Hit_id":
                    hits_vals.append(hit_child.text.split('|')[1]+'_'+hit_child.text.split('|')[2])
                elif hit_child.tag == "Hit_def":
                    hits_vals.append(hit_child.text)
                elif hit_child.tag == "Hit_len":
                    hits_vals.append(hit_child.text)
                elif hit_child.tag == "Hit_hsps":
                    for hsps in hit_child:
                        if hsps.tag == "Hsp":
                            hsps_vals = list()
                            i+=1
                            for hsp_child in hsps:
                                #print('--'+hsp_child.tag+hsp_child.text)
                                if (hsp_child.tag in hsp_cols):
                                    hsps_vals.append(hsp_child.text)
                                hsp_child.clear()
                                while elem.getprevious() is not None:
                                    del elem.getparent()[0]
                                
                            dt_pdb[i] = query_vals+hits_vals+hsps_vals
                            if i % pickle_sz == 0:
                                print(i)
                                if (i==pickle_sz):
                                    mode = 'wb'
                                else:
                                    mode = 'ab'
                                resources.save_pickle(dt_pdb, pdb_filename, mode)
                                dt_pdb = dict()
                        hsps.clear()
                hit_child.clear()
            elem.clear()
    print(i)
    if (mode==""):
        mode = 'wb'
    resources.save_pickle(dt_pdb, pdb_filename, 'wb')
    
    tot_in_sec = time.time() - start_time
    print("--- %s seconds ---" % (tot_in_sec))
    print("Extraction from blast xml file finished.")
    return (pdb_filename)


def count_pdbs_by_eval(dt_pdb):
    ''' Gets the number os different group counts to help on validation. '''
    sum_10e_05 = 0
    sum_001 = 0
    sum_10 = 0
    for key, val in dt_pdb.items():
        if float(val[6]) <= 10e-05:
            sum_10e_05 += 1
        elif (float(val[6]) > 10e-05)&(float(val[6])<=0.01):
            sum_001 += 1
        else:
            sum_10 += 1
    print('# PDBs <= 10e-05: {0}'.format(sum_10e_05))
    print('# PDBs > 10e-05 & <= 0.01: {0}'.format(sum_001))
    print('# PDBs > 0.01 & <= 10: {0}'.format(sum_10))
    

##### FROM dt_pdb #####
def positions2array(df_idr, pdb_path, ids_pos=(0,8,9), cols_df = ['idr_start', 'idr_end'], pickle_sz=10e+05, f_name=''):
    ''' Replicates the start and end positions of each IDR and PDB regions
    to a numpy array.'''
    
    f_name = resources.gen_filename(pdb_path, f_name, "", "pickle")
    
    starts_ends, idx_million = [], []
    len_startsends, i = 0, 0
    IDR_old=''
    mode=''
    start_time = time.time()
    with open(pdb_path, 'rb') as handle:
        while 1:
            try:
                dt_pdb = pickle.load(handle)
                print('loaded: '+str(i))
                for key, val in dt_pdb.items():
                    if (IDR_old) != val[ids_pos[0]]:
                        filtered = df_idr[df_idr['seq_name']==val[ids_pos[0]]]
                        vect = np.array(filtered.loc[:, [cols_df[0], cols_df[1]]])
                    #if len(filtered)>0:
                    rep = np.tile([int(val[ids_pos[1]]), int(val[ids_pos[2]])], (len(filtered),1))
                    pos = np.append(vect, rep, axis=1)
                    starts_ends.append(pos)
                    if ((i>0)&(i%pickle_sz==0)):
                        starts_ends = np.concatenate(starts_ends, axis=0)
                        len_startsends += len(starts_ends)
                        idx_million.append(len_startsends)
                        if (pickle_sz==i):
                            mode = 'wb'
                        else:
                            mode = 'ab'
                        resources.save_pickle(starts_ends, f_name, mode)
                        starts_ends = []
                    IDR_old = val[ids_pos[0]]
                    i+=1
                    if i%10e+04==0:
                        print(i)
            except EOFError:
                break
    print(i)
    starts_ends = np.concatenate(starts_ends, axis=0)
    len_startsends += len(starts_ends)
    idx_million.append(len_startsends)
    if (mode==""):
        mode="wb"
    resources.save_pickle(starts_ends, f_name, mode)
    tot_in_sec = time.time() - start_time
    print("--- %s seconds ---" % (tot_in_sec))
    print("Pairing starts and ends step finished.")
    return idx_million

# un_ids = [0:seq_name, 1:idr_name, 2:hit_id, 3:pdb_name, 4:hsp_id, 5:internalID]
# gen_info = [0:bitscore, 1:hit_len, 2:pdb_start, 3:pdb_end, 4:identity, 5:seq_align_start, 6:seq_align_end, 7:seq_size, 8:internalID, 9:internalID2]
def get_another_data(df_idr, pdb_path, idx_million, ids_pos=(0,1,2,5,6,4,10,11,12,7,14,15), cols_df = ['idr_name', 'idr_size'], pickle_sz=10e+05, dt=(True,True), file_names=''):
    ''' Generates a numpy array of unique IDs and other relevant data and 
    save them to pickle files in disk. '''
    
    assert (sum(dt)==1 or sum(dt)==2), "Please inform a tuple with the values you want to generate (ids and evalues, text data)."
    assert ((sum(dt)==2 and len(file_names)==5) or (sum(dt)==1 and len(file_names)>=4)), "Please provide the file names to be used in the pickles."
    
    sz = idx_million[0]
    mode=""
    if dt[0]:
        evalues = np.empty(sz, dtype=float)
        idr_sizes = np.empty(sz, dtype=int)
        un_ids = np.empty((sz,6), dtype=object)
        gen_info = np.empty((sz,10), dtype=object)
    if dt[1]:
        seqs_blast = np.empty((sz,4), dtype=object)
    i,k,l,m = 0,0,0,0
    IDR_old=''
    start_time = time.time()
    with open(pdb_path, 'rb') as handle:
        while 1:
            try:
                dt_pdb = pickle.load(handle)
                for key, val in dt_pdb.items():
                    if (IDR_old) != val[ids_pos[0]]:
                        filtered = df_idr[df_idr['seq_name']==val[ids_pos[0]]]
                        if dt[0]:
                            ids_idr = filtered.loc[:, ['seq_name', cols_df[0]]].values
                            idr_size = filtered.loc[:, cols_df[1]].values
                            seq_size = filtered.loc[:, 'seq_len'].values
                    
                    # Extracting and merging the IDs to create a unique identifier
                    internal_id = np.arange(l,l+len(filtered)).reshape(len(filtered),1)
                    if dt[0]:     
                        pdb_ids = np.tile([val[ids_pos[1]], val[ids_pos[2]], val[ids_pos[3]]], (len(filtered),1))
                        ids = np.append(ids_idr, pdb_ids, axis=1)
                        ids = np.append(ids, internal_id, axis=1)
                        un_ids[k:k+len(pdb_ids),:] = ids
                        
                        other_pdb = np.tile([val[ids_pos[4]], val[ids_pos[5]], val[ids_pos[6]], val[ids_pos[7]], val[ids_pos[8]], val[ids_pos[13]], val[ids_pos[14]]], (len(filtered),1))
                        internal_id2 = np.repeat(i+1,len(filtered)).reshape(len(filtered),1)
                        #other = np.append(ids_idr, other_pdb, axis=1)
                        other = np.append(other_pdb, seq_size.reshape(len(filtered),1), axis=1) 
                        other = np.append(other, internal_id, axis=1)
                        other = np.append(other, internal_id2, axis=1)
                        gen_info[k:k+len(other_pdb),:] = other
                        
                        # Getting e-values e IDR sizes to filter by homology relevance
                        evalues[k:k+len(other_pdb)] = np.repeat(float(val[ids_pos[9]]),len(filtered))
                        idr_sizes[k:k+len(other_pdb)] = idr_size
                        
                        k=k+len(filtered)
                        l=l+len(filtered)
                        # Saving the values to a pickle to save memory
                        if ((i>0)&(i%pickle_sz==0)):
                            #if (i==(pickle_sz*2)):
                            #    break
                            if (pickle_sz==i):
                                mode = 'wb'
                            else:
                                mode = 'ab'
                            resources.save_pickle(evalues, resources.gen_filename(pdb_path, file_names[0], "", "pickle"), mode)
                            resources.save_pickle(idr_sizes, resources.gen_filename(pdb_path, file_names[1], "", "pickle"), mode)
                            resources.save_pickle(un_ids, resources.gen_filename(pdb_path, file_names[2], "", "pickle"), mode)
                            resources.save_pickle(gen_info, resources.gen_filename(pdb_path, file_names[3], "", "pickle"), mode)
                            if len(idx_million) > i/pickle_sz:
                                sz = idx_million[int(i/pickle_sz)]-idx_million[int((i/pickle_sz)-1)]
                                evalues = np.empty(sz, dtype=float)
                                idr_sizes = np.empty(sz, dtype=int)
                                un_ids = np.empty((sz,6), dtype=object)
                                gen_info = np.empty((sz,10), dtype=object)
                                k=0
                            
                    # Extracting the blast overlaps and midline
                    if dt[1]:
                        blast_pdb = np.tile([val[ids_pos[10]], val[ids_pos[11]], val[ids_pos[12]]], (len(filtered),1))
                        blast = np.append(internal_id, blast_pdb, axis=1) 
                        seqs_blast[m:m+len(blast_pdb),:] = blast
                        m=m+len(filtered)
                        
                        if ((i>0)&(i%pickle_sz==0)):
                            resources.save_pickle(seqs_blast, resources.gen_filename(pdb_path, file_names[4], "", "pickle"), mode)
                            if len(idx_million) > i/pickle_sz:
                                sz = idx_million[int(i/pickle_sz)]-idx_million[int((i/pickle_sz)-1)]
                                seqs_blast = np.empty((sz,4), dtype=object)
                                m=0
                    
                    IDR_old = val[0]
                    i+=1
                    if i%10e+04==0:
                        print(i)
            except EOFError:
                break
    print(i)
    if (mode==""):
        mode="wb"
    if dt[0]:
        resources.save_pickle(evalues, resources.gen_filename(pdb_path, file_names[0], "", "pickle"), mode)
        resources.save_pickle(idr_sizes, resources.gen_filename(pdb_path, file_names[1], "", "pickle"), mode)
        resources.save_pickle(un_ids, resources.gen_filename(pdb_path, file_names[2], "", "pickle"), mode)
        resources.save_pickle(gen_info, resources.gen_filename(pdb_path, file_names[3], "", "pickle"), mode)
    if dt[1]:
        resources.save_pickle(seqs_blast, resources.gen_filename(pdb_path, file_names[4], "", "pickle"), mode)
    tot_in_sec = time.time() - start_time
    print("--- %s seconds ---" % (tot_in_sec))
    print("Extraction from blast pickle finished.")
    

##### PART 2, CROSSING IDRS AND PDBs AND EXTRACTING ALL VIABLE COMBINATIONS #####
    
##### AUXILIARY FUNCTIONS get_all #####
def define_overlaps(pos):
    ''' Gets the start and end to the regions with overlaps and set 0 to the ones 
    without them. '''
    #pos = np.array([[22,55,139,237], [40,73,64,96],[62,87,12,22], [64,96,73,82], [73,82,64,96], [64,96,40,73]])
    idx_mat = np.argsort(pos, axis=1)
    r1 = np.transpose(np.tile((pos[:,1]-pos[:,2]>0) & (pos[:,3]-pos[:,0]>0), (4,1)))
    ord_mat = np.sort(pos, axis=1)
    nrow = len(ord_mat)
    r2 = np.concatenate((np.zeros((nrow,1)), np.ones((nrow,2)), np.zeros((nrow,1))), axis=1)
    #r2_2 = np.transpose(np.tile((pos[:,0]>pos[:,2]) & (pos[:,3]>pos[:,1]), (4,1)))
    #r2 = np.logical_xor(r2_1, r2_2)
    res = np.multiply(np.multiply(ord_mat, r1), r2)
    res = np.sort(res, axis=1)[:,2:]
    sz_overlaps = res[:,1]-res[:,0]
    sz_overlaps[sz_overlaps>0] = sz_overlaps[sz_overlaps>0]+1
    return res, sz_overlaps


def filter_kept_idrs(idrs_removed, un_ids, starts_ends, idr_sizes, evalues, gen_info):
    ''' Removing IDRs with overlaps already assigned to a lower homology group. 
    This action will both save memory and make sure no data be wrongly overlaped. '''
    if len(idrs_removed)>0:
        _, pdb_removed0, idx_all_un = np.unique(un_ids[:,1], return_index=True, return_inverse=True)
        _, pdb_removed1, _ = np.intersect1d(un_ids[:,1], idrs_removed, return_indices=True)
        idx_rep_removed = np.isin(pdb_removed0[idx_all_un], pdb_removed1)
        un_ids_n = np.delete(un_ids, idx_rep_removed, 0)
        starts_ends_n = np.delete(starts_ends, idx_rep_removed, 0)
        idr_sizes_n = np.delete(idr_sizes,  idx_rep_removed, 0)
        evalues_n = np.delete(evalues, idx_rep_removed, 0)
        gen_info_n = np.delete(gen_info, idx_rep_removed, 0)
    return un_ids_n, starts_ends_n, idr_sizes_n, evalues_n, gen_info_n
    

def filter_no_pdb(df_idrs_dt, un_ids):
    ''' Extracts a list of the IDRs with no overlap with PDB. 
    Returns:
        a_idrs_nopdb = list of unique IDR IDs with no alignment with PDB;
        a_prots_nopdb = list of uniprot IDs with no overlap with PDB.
    '''
    a_idrs_un = np.array(df_idrs_dt['seq_name'].unique().tolist())
    a_pdbs_un = np.unique(un_ids[:,0])
    a_seqs_notpdb = np.setdiff1d(a_idrs_un, a_pdbs_un)
    a_idrs_nopdb = df_idrs_dt[df_idrs_dt['seq_name'].isin(a_seqs_notpdb)]
    a_prots_nopdb = a_idrs_nopdb['seq_name'].unique().tolist()
    return a_idrs_nopdb, a_prots_nopdb

#Control Q9BUP0 - short with long
#Control Q92793 - short with some long
def filter_no_overlap(df_idrs_dt,un_ids,pos_over,prefix):
    ''' Filter the cases where homologous sequences were found, but the IDR
    region did not overlap. '''
    # Unique IDR IDs with counts (they were doubled by the PDB occurrences)
    total_counts = np.unique(un_ids[:, 1], return_counts=True)
    # Unique IDR IDs with counts for the ones without overlaps
    total_nooverlap = np.unique(un_ids[pos_over[:,0]==0, 1], return_counts=True)
    # Gets the intersection between both groups, returning the first occurrence index
    inters = np.intersect1d(total_counts[0],total_nooverlap[0], assume_unique=True, return_indices=True)
    # Filters the IDR names and counts by first occurrence index
    idrs_counts_filt = (total_counts[0][inters[1]], total_counts[1][inters[1]])
    # gets the IDR names where no PDB had overlaps (because some of the PDB may have overlaps and some may not)
    ids_nooverlap = total_nooverlap[0][total_nooverlap[1]==idrs_counts_filt[1]]
    # Filters the Dataframe to keep just the IDRs with no overlaps at all
    idrs_nooverlap = df_idrs_dt[df_idrs_dt[prefix+'_name'].isin(ids_nooverlap)]
    #prots_nooverlap = idrs_nooverlap['seq_name'].unique().tolist()
    return idrs_nooverlap#, prots_nooverlap


def get_overlap_size(sz_over, idr_sizes):
    ''' Calculates the fraction of the region overlaped with a PDB sequence. '''
    over_perc = sz_over/idr_sizes
    return over_perc


def extract_short_overlaps(sz_over, over_perc, cutoff=.5, min_sz=30):
    ''' Get the short overlaps bool list and the not long ones, to ensure no
    short are selected if there are long ones in the same IDR. Returns a boolean
    array to use in the selection of the short and long ones. '''
    check_short = ((over_perc<cutoff)&(sz_over<min_sz)&(over_perc>0))
    check_long = (over_perc>0)&((sz_over>=min_sz)|(over_perc>=cutoff))
    return check_short, check_long


def filter_short_overlaps(df_idrs_dt,un_ids,short_bool,long_bool,prefix):
    ''' Filter the cases where homologous sequences were found, but the IDR
    region had a short overlap. '''
    pdbs_long = np.unique(un_ids[long_bool, 1])
    short_over = np.unique(un_ids[short_bool, 1])
    short_over = short_over[~np.isin(short_over, pdbs_long)]
    idrs_short_over = df_idrs_dt[df_idrs_dt[prefix+'_name'].isin(short_over)]
    return idrs_short_over


def filter_overlaps(df_idrs_dt,un_ids,pos_over,starts_ends,evalues,gen_info,short_bool,long_bool,over_perc,sz_over,cols,prefix):
    ''' Gets all the base data and filter the long and short overlaps. 
    Outputs:
        df_idrs_dt_over = The subset of the IDR dataframe with actual overlaps;
        pdb_over = A dataframe with all data related to the PDB overlap;
        pdb_short = A dataframe with all data related to short overlaps only (for validation purposes).
    '''
    p_ids_all = un_ids[long_bool, :-1]
    p_posit_over = pos_over[long_bool, :]
    p_posit_orig = starts_ends[long_bool, :]
    p_evals = evalues[long_bool]
    p_overl_perc = over_perc[long_bool]
    p_sz_over = sz_over[long_bool].reshape(len(p_overl_perc),1)
    p_gen_info = gen_info[long_bool, :]
    pdb_over = pd.DataFrame(np.concatenate([p_ids_all,p_gen_info,p_posit_over,p_posit_orig,p_evals.reshape(len(p_evals),1), p_overl_perc.reshape(len(p_overl_perc),1), p_sz_over], axis=1), columns=cols)

    p_ids_all = un_ids[short_bool, :-1]
    p_posit_over = pos_over[short_bool, :]
    p_posit_orig = starts_ends[short_bool, :]
    p_evals = evalues[short_bool]
    p_overl_perc = over_perc[short_bool]
    p_sz_over = sz_over[short_bool].reshape(len(p_overl_perc),1)
    p_gen_info = gen_info[short_bool, :]
    pdb_short = pd.DataFrame(np.concatenate([p_ids_all,p_gen_info,p_posit_over,p_posit_orig,p_evals.reshape(len(p_evals),1), p_overl_perc.reshape(len(p_overl_perc),1), p_sz_over], axis=1), columns=cols)
    
    idrs_over = np.unique(un_ids[long_bool, 1])
    df_idrs_dt_over = df_idrs_dt[df_idrs_dt[prefix+'_name'].isin(idrs_over)]
    return df_idrs_dt_over, pdb_over, pdb_short


def filter_evalue(evalues, evalue_min=None, evalue_max=10e-05):
    ''' Cut the e-value array to keep just the rows in the e-value range. '''
    if evalue_min==None:
        eval_cut = evalues<=evalue_max
    else:
        eval_cut = (evalues>evalue_min)&(evalues<=evalue_max)
    return eval_cut
    

def filter_pdb_array(arr_data, filtered_eval):
    ''' Filter the data extracted from the pickle file to keep just the rows
    in the selected e-value range.
    '''
    if arr_data.ndim==1:
        filtered_arr = arr_data[filtered_eval]
    else:
        filtered_arr = arr_data[filtered_eval,:]
    return filtered_arr


def merge_noIDRs_pdbs(df_noPDB, df_idrs_dt, pdbs_overlap, min_eval, max_eval, prefix):
    ''' Adds the PDB overlaps data to the IDR dataframe. '''
    
    if max_eval==0.0001:
        smax_eval = 'b_10e-05'
        smax_eval_name = '10e-05'
    elif max_eval==0.01:
        smax_eval = 'c_0.01'
        smax_eval_name = '0.01'
    else:
        smax_eval = 'd_10'
        smax_eval_name = '10'
        
    if min_eval==0.0001:
        smin_eval = 'b_10e-05'
    elif min_eval==0.01:
        smin_eval = 'c_0.01'
    elif min_eval==None:
        smin_eval = 'a_Overlap'
    else:
        smin_eval = 'd_10'
    
    cols2add = [prefix+'_name', 'pdb_name', 'pdb_evalue', 'over_perc', 'over_sz', 'bitscore', 'internal_id', 'internal_id2', 'seq_start', 'seq_end', 'seq_align_start', 'seq_align_end', 'pdb_start', 'pdb_end', 'hit_id', 'hsp_id', 'identity']
    df_all_pdbs = pdbs_overlap.loc[:, cols2add]
    df_all_pdbs.loc[:, 'pdb_name'] = df_all_pdbs.loc[:, 'pdb_name'].astype(str)
    for cols in cols2add[2:]:
        df_all_pdbs.loc[:, cols] = df_all_pdbs.loc[:, cols].astype(float)
    
    # End up with more than a thousand alignments to some proteins
    #df_all_pdbs = df_all_pdbs.sort_values(by=[prefix+'_name', 'over_perc','pdb_evalue', 'bitscore'], ascending=[True, False, True, False]).drop_duplicates(subset=[prefix+'_name'])
    df_all_pdbs = df_all_pdbs.sort_values(by=[prefix+'_name', 'over_perc','pdb_evalue', 'bitscore'], ascending=[True, False, True, False])
    # When I filtered to 20, in some cases all alignments where synthetic.
    # Evaluating without the filter if this changes.
    #df_all_pdbs["val"] = 1
    #df_all_pdbs["num_count"] = df_all_pdbs.groupby(prefix+'_name')['val'].cumsum()
    #df_all_pdbs = df_all_pdbs.loc[df_all_pdbs["num_count"]<=20, cols2add]
    
    df_all_pdbs = df_all_pdbs.loc[~df_all_pdbs[prefix+'_name'].isin(df_noPDB[prefix+'_name'].unique())]
    df_all_pdbs['removed_on_s'] = smax_eval
    df_all_pdbs['removed_on'] = smax_eval_name
    df_all_pdbs['kept_on'] = smin_eval
    df_all_pdbs['kept_on'] = df_all_pdbs['kept_on'].astype(str)
    
    #new_df_idrs_dt = df_idrs_dt.copy()
    sel = df_idrs_dt[prefix+'_name'].isin(df_all_pdbs[prefix+'_name'].unique())
    if 'pdb_name' in df_idrs_dt:
        sel_cols = cols2add[1:]+['removed_on_s', 'removed_on', 'kept_on']
        add_df_idrs = df_idrs_dt.loc[sel, df_idrs_dt.columns.difference(sel_cols)]
        add_df_idrs = pd.merge(add_df_idrs, df_all_pdbs, how='left', on=prefix+'_name')
        df_idrs_dt = df_idrs_dt.drop(df_idrs_dt[sel].index)
        df_idrs_dt = pd.concat([df_idrs_dt, add_df_idrs], ignore_index=True)
        #for cols in cols2add[1:]:
        #    df_idrs_dt.loc[sel, cols] = list(df_all_pdbs.loc[:, cols])
    else:
        df_idrs_dt = pd.merge(df_idrs_dt, df_all_pdbs, how='left', on=prefix+'_name')
    df_idrs_dt = df_idrs_dt.sort_values(by=['seq_name', prefix+'_num'])
    return df_idrs_dt


def group_min(groups, data):
    ''' Extracts the minimum value by group from unsorted values. '''
    order = np.lexsort((data, groups))
    groups = groups[order]
    data = data[order]
    index = np.empty(len(groups), 'bool')
    index[0] = True
    index[1:] = groups[1:] != groups[:-1]
    return data[index]


def group_max(groups, data):
    ''' Extracts the maximum value by group of unsorted values. '''
    order = np.lexsort((data, groups))
    groups = groups[order] #this is only needed if groups is unsorted
    data = data[order]
    index = np.empty(len(groups), 'bool')
    index[-1] = True
    index[:-1] = groups[1:] != groups[:-1]
    return data[index]


def eval_no_overlaps(over_perc, gen_info, starts_ends, un_ids, evalues, prefix):
    ''' Generates some statistical analysis based on the cases without 
    PDB overlaps to support further analysis. '''
    bool_nopdb = over_perc==0
    seq_lens = gen_info[bool_nopdb, 5].astype(int)
    evalues_pdb = evalues[bool_nopdb]
    start_pdb = starts_ends[bool_nopdb, 2]
    end_pdb = starts_ends[bool_nopdb, 3]
    seq_ids = un_ids[bool_nopdb, 1]
    seq_un, seq_un_idx, rep_idx = np.unique(seq_ids, return_index=True, return_inverse=True)
    seq_lens = seq_lens[seq_un_idx]
    start_min = group_min(rep_idx, start_pdb)
    start_rel_min = start_min/seq_lens
    end_min = group_min(rep_idx, end_pdb)
    end_rel_min = end_min/seq_lens
    start_max = group_max(rep_idx, start_pdb)
    start_rel_max = start_min/seq_lens
    end_max = group_max(rep_idx, end_pdb)
    end_rel_max = end_max/seq_lens
    evalue_min = group_min(rep_idx, evalues_pdb)
    evalue_max = group_max(rep_idx, evalues_pdb)
    cols = [prefix+"_name", "nopdb_start_min", "nopdb_end_min", "nopdb_rel_start_min", 
            "nopdb_rel_end_min", "nopdb_start_max", "nopdb_end_max", 
            "nopdb_rel_start_max", "nopdb_rel_end_max", "nopbdb_evalue_min",
            "nopbdb_evalue_max"]
    df_info_nopdb = pd.DataFrame(np.hstack((seq_un[:,None], start_min[:,None], 
                                            end_min[:,None], start_rel_min[:,None], 
                                            end_rel_min[:,None], start_max[:,None], 
                                            end_max[:,None], start_rel_max[:,None],
                                            end_rel_max[:,None], evalue_min[:,None], 
                                            evalue_max[:,None])), columns=cols)
    return df_info_nopdb
    

def report_sizes(idrs_nopdb, idrs_nooverlap, idrs_shortoverlap, df_idrs_dt, prots_nopdb, idrs_all_nopdb, evals):
    ''' A simple report just to give us an overview of how many overlaps actually
    happened. The idea is that the information about the cases without homologous
    at all can be usefull. The counts are on each iteration, however, not for 
    all e-values at once. '''
    print()
    print('Evalues equal/smaller than {0}: '.format(evals))
    print('Total IDRs: {0}'.format(len(df_idrs_dt)))
    prots = len(df_idrs_dt['seq_name'].unique())
    print('Total of Proteins with IDRs: {0}'.format(prots))
    perc = round((len(prots_nopdb)/prots)*100,2)
    print('Proteins with NO Homologous AT ALL: {0} ({1}%)'.format(len(prots_nopdb), perc))
    perc = round((len(idrs_nopdb)/len(df_idrs_dt))*100, 2)
    print('IDRs without Homologous in the PDB: {0} ({1}%)'.format(len(idrs_nopdb), perc))
    perc = round((len(idrs_nooverlap)/len(df_idrs_dt))*100, 2)
    print('IDRs without Overlaps in the PDB: {0} ({1}%)'.format(len(idrs_nooverlap), perc))
    perc = round((len(idrs_shortoverlap)/len(df_idrs_dt))*100, 2)
    print('IDRs without short Overlaps (up to 10 AAs or 0.05%): {0} ({1}%)'.format(len(idrs_shortoverlap), perc))
    perc = round((len(idrs_all_nopdb)/len(df_idrs_dt))*100, 2)
    print('Total IDRs NO PDB: {0} ({1}%)'.format(len(idrs_all_nopdb), perc))
    prots_no = len(idrs_all_nopdb['seq_name'].unique())
    perc = round((prots_no/prots)*100,2)
    print('Total Proteins NO PDB: {0} ({1}%)'.format(prots_no, perc))
    print()

##### FUNC get_all #####
def generate_all_dfs(basis_path, df_idr, min_eval, max_eval, prefix="idr", cutoff=(.5,30), save=False, file_names=''):
    ''' Extract the information of the pair IDR/PDB and add to a global dataframe
    
    Inputs:
        basis_path = path where the pickle files are saved;
        df_idr = dataframe with all the available IDR data. New columns will be added
                to this dataframe, returning a more complete file with overlaps info;
        min_eval = minimum e-value for the round. The main goal is to make sure 
                    that the IDRs with lower e-value get assigned again with an
                    alignment to a PDB sequence with a higher e-value;
        max_eval = maximum e-value for the round. The main goal is to make sure 
                    that the IDRs with lower e-value get assigned again with an
                    alignment to a PDB sequence with a higher e-value;
        prefix = The function was designed to run over different kinds of input,
                    as long as the pickle files got the proper starts-end positions
                    informed, you can use any partition of the sequence, not only
                    IDRs. The prefix define the names of the resulting columns,
                    e.g. idr, poly, polyxy...
        cutoff = Define the limits of the overlaps to give some flexibility for
                    cases not 100% aligned. The default accepts 50% of overlap 
                    up to a minimum of 30 amino acids. This criteria was defined
                    based on the smallest size of an IDR in mobidb, 20 amino acids;
        save = Additional validation files should be saved to disk?
        file_names = List of names for the pickle files. An exact order must be
                    respected to make sure the function works properly. Defined
                    as a parameter because this files will be used several times
                    in the future.
                    
    Outputs:
        df_idr = The original IDR file with the extra overlapping columns. As we
                need to perform an extra step to filter 1 PDB per IDR, this file
                will still have several PDBs per IDR. We realized that somethimes 
                the alignment would still bring IDRs overlaping complete missing 
                regions with good e-value, so we moved the selection criteria 
                for after we map the SS data and remove synthetic alignments as well;
        idrs_all_nopdb = Original df_idr file with the cases with no overlap with
                PDB or with overlaps shorter than cutoff. File generated for 
                validation purposes;
        pdbs_short = List of PDB alignments with short overlap with IDRs. This
                cases are not available anywhere else after the filtering. File
                generated for validation purposes;
        idrs_overlap (save): List of IDRs with overlaps independently if they
                are short or long;
        df_info_nopdb (save): Statistical data related to the cases with no 
                overlap to PDB (validation purposes).
    '''
    # Loading pickles from disk
    starts_ends = resources.open_pickle(file_names[0], basis_path)
    evalues = resources.open_pickle(file_names[1], basis_path)
    idr_sizes = resources.open_pickle(file_names[2], basis_path)
    un_ids = resources.open_pickle(file_names[3], basis_path)
    gen_info = resources.open_pickle(file_names[4], basis_path)
    
    # Removing the IDRs with overlaps in the last iteration our strategy appends
    # the cases with overlaps in the iteration to the previous one
    if ('removed_on' in df_idr) and (not save):
        df_idr_filtered = df_idr.loc[df_idr['removed_on'].isna()]
        idrs_removed = df_idr.loc[~(df_idr['removed_on'].isna()), prefix+'_name'].unique()
        un_ids, starts_ends, idr_sizes, evalues, gen_info = filter_kept_idrs(idrs_removed, un_ids, starts_ends, idr_sizes, evalues, gen_info)
    else:
        df_idr_filtered = df_idr.copy()
    
    # Filters the data extracted from pickles based on the e-value range
    pdb_qual = filter_evalue(evalues, min_eval, max_eval)
    starts_ends = filter_pdb_array(starts_ends, pdb_qual)
    un_ids = filter_pdb_array(un_ids, pdb_qual)
    idr_sizes = filter_pdb_array(idr_sizes, pdb_qual)
    evalues = filter_pdb_array(evalues, pdb_qual)
    gen_info = filter_pdb_array(gen_info, pdb_qual)
    # Generates overlap information
    pos_over, sz_over = define_overlaps(starts_ends)
    over_perc = get_overlap_size(sz_over, idr_sizes)
    short_bool, long_bool = extract_short_overlaps(sz_over, over_perc, cutoff[0], cutoff[1])
    #del locals()['starts_ends']
    
    # Distinguish IDRs without homologs at all and the ones with short overlap
    # from the cases that simply don't overlap
    # These cases will be concatenated in one dataframe as not overlaped with 
    # a column to use as a class.
    idrs_nopdb, prots_nopdb = filter_no_pdb(df_idr_filtered,un_ids)
    idrs_nopdb = idrs_nopdb.assign(pdb_type='No Homologs')
    idrs_nooverlap = filter_no_overlap(df_idr_filtered,un_ids,pos_over,prefix)
    idrs_nooverlap = idrs_nooverlap.assign(pdb_type='No Overlaps')
    # Filter IDRs with short overlap - up to 50% coverage or 10 AA
    idrs_shortoverlap = filter_short_overlaps(df_idr_filtered,un_ids,short_bool,long_bool, prefix)
    idrs_shortoverlap = idrs_shortoverlap.assign(pdb_type='Short Overlaps')
    idrs_all_nopdb = pd.concat([idrs_nopdb, idrs_nooverlap, idrs_shortoverlap], ignore_index=True)
    
    # Filtering IDRs with long overlaps (for validation purposes)
    cols=['seq_name', prefix+'_name', 'hit_id', 'pdb_name', 'hsp_id', 'bitscore', 
          'hit_len', 'pdb_start', 'pdb_end', 'identity', 'seq_align_start', 
          'seq_align_end', 'seq_len', 'internal_id', 'internal_id2', 'over_start', 'over_end', 
          prefix+'_start', prefix+'_end', 'seq_start', 'seq_end', 'pdb_evalue', 'over_perc', 'over_sz']
    idrs_overlap, pdbs_overlap, pdbs_short = filter_overlaps(df_idr_filtered,un_ids,pos_over,starts_ends,evalues,gen_info,short_bool,long_bool,over_perc,sz_over,cols,prefix)
    report_sizes(idrs_nopdb, idrs_nooverlap, idrs_shortoverlap, df_idr, prots_nopdb, idrs_all_nopdb, max(evalues))
    df_info_nopdb = eval_no_overlaps(over_perc, gen_info, starts_ends, un_ids, evalues, prefix)
    if save:
        name = re.sub(r'\.', '-', str(max_eval))
        idrs_all_nopdb.to_csv(basis_path+'nopdb_all_'+prefix+'_'+name+'.csv', index=False)
        idrs_overlap.to_csv(basis_path+'df_over_'+prefix+'_'+name+'.csv', index=False)
        pdbs_short.to_csv(basis_path+'df_shortpdb_'+prefix+'_'+name+'.csv', index=False)
        df_info_nopdb.to_csv(basis_path+'nopdb_info_'+prefix+'_'+name+'.csv', index=False)
    else:
        df_idr = merge_noIDRs_pdbs(idrs_all_nopdb, df_idr, pdbs_overlap, min_eval, max_eval, prefix)
        
    return df_idr, idrs_all_nopdb, pdbs_short

##### PART 3, GETTING THE 2D STRUCTURE ANNOTATIONS AND FILTERING THE BEST IDR/PDB PAIRS #####

##### Auxiliary functions

def mark_no_homologs(basis_path, df_idr, save_name):
    ''' Re-check to make sure the class of each filtering groups is properly 
    assigned because this variable will be used to filter the significant 
    alignemnts later'''
    df_idr.loc[df_idr['removed_on_s'].isna(), 'removed_on_s'] = 'e_no_homology'
    df_idr.loc[df_idr['removed_on'].isna(), 'removed_on'] = 'No Homology'
    df_idr.to_csv(basis_path+save_name, index=False)
    return(df_idr)


def append_seq_details(seqs_path, df_idr_details):
    ''' Gets additional sequence data that will be used on the 2D structure annotations.
    This data could have being added before, but we just realized later on in the process
    the need for it. Decided to add it when needed.'''
    seq_dict = resources.load_seqs(seqs_path)
    df_idr_details['seq_desc'] = df_idr_details["seq_name"].apply(lambda x: seq_dict[x][1])
    df_idr_details['seq_aa'] = df_idr_details["seq_name"].apply(lambda x: seq_dict[x][0])
    return df_idr_details


def extract_pdb_sizes(ss_dict):
    ''' Calculates the size of the complete PDB structure using just the PDB ID
    and adding the sizes of individual chains. '''
    dct_sizes = defaultdict(list)
    for k, v in ss_dict.items():
        dct_sizes[k.split("_")[0]].append(len(v[0]))
    return pd.DataFrame.from_dict({key: sum(dct_sizes[key]) for key in dct_sizes}, columns=["pdb_all_size"], orient="index").reset_index().rename(columns={"index":"pdb_id"})


def append_pdb_details(pdb_mask_path, df_idr_details):
    ''' Gets additional PDB data that will be used on the 2D structure annotations.'''
    ss_dict = resources.load_seqs(pdb_mask_path, "|", 'pdb')
    pdb_sizes = extract_pdb_sizes(ss_dict)
    df_idr_details["pdb_id"] = df_idr_details.pdb_name.str.split("_", n=1, expand=True)[0]
    df_idr_details = pd.merge(df_idr_details, pdb_sizes, how="left", on="pdb_id")
    df_idr_details['pdb_desc'] = df_idr_details["pdb_name"].apply(lambda x: ss_dict[x][1] if(str(x) != 'nan') else x)
    df_idr_details['pdb_aa'] = df_idr_details["pdb_name"].apply(lambda x: ss_dict[x][0] if(str(x) != 'nan') else x)
    df_idr_details['pdb_desc'] = df_idr_details['pdb_desc'].astype(str)
    df_idr_details['pdb_size'] = df_idr_details["pdb_name"].apply(lambda x: len(ss_dict[x][0]) if(str(x) != 'nan') else x)
    df_idr_details['over_perc_seq'] = round((df_idr_details["seq_end"]-df_idr_details["seq_start"])/df_idr_details["seq_len"], 2)
    df_idr_details['over_perc_pdb'] = round((df_idr_details["pdb_end"]-df_idr_details["pdb_start"])/df_idr_details["pdb_size"], 2)
    df_idr_details['over_dir'] = df_idr_details.apply(lambda x: -1 if x["seq_start"] > x[PREF+"_start"] else 1, axis=1)
    return df_idr_details


def append_ss_details(ssSeq_path, df_idr_details):
    ''' Extract secondary structure data (SS) from the file in fasta format 
    composed by 1 entry per sequence and 1 entry for the SS data. '''
    ss_data, _, _, _ = resources.extract_ss(ssSeq_path)
    df_idr_details['ss_seq'] = df_idr_details["pdb_name"].apply(lambda x: ss_data[x] if(str(x) != 'nan') else x)
    return df_idr_details




def merge_idr_pdb(idrs_path, pdb_path, file_names, cutoff):
    ''' Parts 1 and 2 - Main function IDR-PDB, no filtering. '''
    # This process can take several minutes to run.
    
    basis_path = resources.get_dir(pdb_path)+"/"+resources.get_filename(pdb_path)+"_"
    
    # Not sure if the separators can be adjusted in the blast file
    sep1 = "|"
    sep2 = "|"
    # Dumping numpy objects to disk to save memory
    pickle_sz=10e+05
    
    # Starting with the extraction of each blast pair from the xml file
    blast_pickle = extract_blast_pairs(pdb_path, pickle_sz, sep1, sep2)
    df_idr = pd.read_csv(idrs_path)
    
    # Replicating the PDBs and IDRs to extract Positions, E-value and unique IDSs an all other relevant alignment data
    idx_million = positions2array(df_idr, blast_pickle, (0,8,9), ['idr_start', 'idr_end'], pickle_sz, 'starts_ends')
    get_another_data(df_idr, blast_pickle, idx_million, (0,1,2,5,6,4,10,11,12,7,14,15,16,8,9), 
                                   ['idr_name', 'idr_size'], pickle_sz, (1,1), file_names)
    
    if len(file_names)==5:
        file_names.insert(0, 'starts_ends')
    
    # We started running for less significant e-values too, but data proved not relevant
    # So we dropped the other e-values: 0.01 (maybe significant) and 10 (not significant).
    # Check the function to know details about the other outputs. They were created for validation purposes
    df_idr_new, _, _ = generate_all_dfs(basis_path, df_idr, None, 10e-05, "idr", cutoff, False, file_names)
    # Run the same process passing df_idr_new as 2nd parameter to add less significant alignments.
    # If they are not relevant, you may consider add the significance paramenter to local blast 
    # to reduce the number of alignments and improve performance.
    df_idr_new = mark_no_homologs(basis_path, df_idr_new, 'EVAL_data_noSeqs_idr.csv')
    

def run_ss_annotation(idrs_path, pdb_mask_path, ssSeq_path):
    ''' Parts 3 - Main function IDR-PDB, getting the 2D structure and selecting 
    best candidate. '''
    
    basis_path = resources.get_dir(pdb_mask_path)+"/"+resources.get_filename(pdb_mask_path)+"_"
    
    if len(file_names)==5:
        file_names.insert(0, 'starts_ends')
        
    df_idr_new = pd.read_csv(idrs_path)
    # ===================== Get PDB overlaps proportions =====================
    df_idr_details = append_pdb_details(pdb_mask_path, df_idr_new)
    df_idr_details = append_ss_details(ssSeq_path, df_idr_details)
    # ========================== Add seqs info ===============================
    # Need to add sequence details first to validate the real sequence start/end
    # position of the alignment later
    df_idr_details = append_seq_details(seqs_path, df_idr_details)
    df_idr_details, df_2D_idr, pdb_data = get_dssp_idr(dssp_path, basis_path+file_names[5], df_idr_details)
    
