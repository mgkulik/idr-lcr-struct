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
import pdbDssp

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

def apply_pdb_data(df_idr_new, pdb_det_path):
    ''' We decided not to use synthetic PDB structures in the process. There are
    several definitions of what synthetic means, so we are removing the cases
    where the organism of the protein is synthetic, based on the pdb/cif file. 
    We also add the experiment resolution info to the process. '''
    df_pdb_details = pd.read_csv(pdb_det_path)
    # Extracting just the PDB id (without chain) because all the structure is annotated as synthetic
    df_idr_new["pdb_id"] = df_idr_new.pdb_name.str.split("_", n=1, expand=True)[0]
    df_idr_new = pd.merge(df_idr_new, df_pdb_details, how="left", on="pdb_id")
    # The merging process added nan in the column and was affecting the filtering process
    df_idr_new = df_idr_new.loc[df_idr_new["synthetic"]!=True, :]
    df_idr_new.drop('synthetic', axis=1, inplace=True)
    return df_idr_new

#### PART 3, GETTING THE 2D STRUCTURE ANNOTATIONS AND FILTERING THE BEST IDR/PDB PAIRS ####

#### AUXILIARY FUNCTIONS ####

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
    df_idr_details = pd.merge(df_idr_details, pdb_sizes, how="left", on="pdb_id")
    df_idr_details['pdb_desc'] = df_idr_details["pdb_name"].apply(lambda x: ss_dict[x][1] if(str(x) != 'nan') else x)
    df_idr_details['pdb_aa'] = df_idr_details["pdb_name"].apply(lambda x: ss_dict[x][0] if(str(x) != 'nan') else x)
    df_idr_details['pdb_desc'] = df_idr_details['pdb_desc'].astype(str)
    df_idr_details['pdb_size'] = df_idr_details["pdb_name"].apply(lambda x: len(ss_dict[x][0]) if(str(x) != 'nan') else x)
    df_idr_details['over_perc_seq'] = round((df_idr_details["seq_end"]-df_idr_details["seq_start"])/df_idr_details["seq_len"], 2)
    df_idr_details['over_perc_pdb'] = round((df_idr_details["pdb_end"]-df_idr_details["pdb_start"])/df_idr_details["pdb_size"], 2)
    df_idr_details['over_dir'] = df_idr_details.apply(lambda x: -1 if x["seq_start"] > x["idr_start"] else 1, axis=1)
    return df_idr_details


def append_ss_details(ssSeq_path, df_idr_details):
    ''' Extract secondary structure data (SS) from the file in fasta format 
    composed by 1 entry per sequence and 1 entry for the SS data. '''
    ss_data, _, _, _ = resources.extract_ss(ssSeq_path)
    df_idr_details['ss_seq'] = df_idr_details["pdb_name"].apply(lambda x: ss_data[x] if(str(x) != 'nan') else x)
    return df_idr_details


#### SET OF FUNCTIONS TO GET 2D STRUCTURE SPECIFIC REGION AND SURROUNDINGS ####


def insert_gaps_2D(seq, struc2D, tp=1, ch1='-', ch2='-'):
    ''' Manages the gaps from the alignment depending on procedure should be executed:   
        tp=1 replace the value in the position
        tp=2 add an extra value to the position. '''
    pos = [i for i, letter in enumerate(seq) if letter == ch1]
    for c in pos:
        if (tp==1):
            struc2D = struc2D[0:c]+ch2+struc2D[c:]
        else:
            struc2D = struc2D[0:c]+ch2+struc2D[c+1:]
    return(struc2D)


def transform_midlines(midline, align_seq, align_pdb):
    ''' BlastP midline is too confusing. Changing it to ClustalOmega patterns. '''
    midline_new = insert_gaps_2D(align_seq, midline, 2, "X", "-")
    midline_new = insert_gaps_2D(align_seq, midline_new, 2, "-", "-")
    midline_new = insert_gaps_2D(align_pdb, midline_new, 2, "X", "-")
    midline_new = insert_gaps_2D(align_pdb, midline_new, 2, "-", "-")
    midline_mask = re.sub(r'\+', ":", midline_new)
    midline_mask = re.sub(r'\w', "*", midline_mask)
    midline_mask = re.sub(r'\s', ".", midline_mask)
    midline_mask = re.sub(r'\-', " ", midline_mask)
    return(midline_mask)


def get_blast_seqs(df_basis, blast_over_path):
    ''' Get extra information about regions aligned and midline generated by
    blastP. Loaded only when really necessary.'''
    reduced=None
    accum_ant=-1
    #handle = open(blast_over_path, 'rb')
    #blast_over = pickle.load(handle)
    with open(blast_over_path, 'rb') as handle:
        while 1:
            try:
                big_list_ids = df_basis.loc[~df_basis['internal_id'].isna(),'internal_id'].astype(int).values
                blast_over = pickle.load(handle)
                big_list_ids = big_list_ids[(big_list_ids>accum_ant)&(big_list_ids<len(blast_over)+accum_ant)]
                big_list_ids = big_list_ids-(accum_ant+1)
                reduced_part = pd.DataFrame(blast_over[big_list_ids,:], columns = ['internal_id', 'align_seq', 'align_pdb', 'midline'])
                if reduced is None:
                    reduced = reduced_part
                else:
                    if len(reduced_part)>0:
                        reduced = pd.concat([reduced, reduced_part], ignore_index=True)
                # As I started the variable with -1, I'm deducing it and when 
                # comparing the sizes I'm using bellow and not bellow-equal
                accum_ant += len(blast_over)
            except EOFError:
                break
    reduced['internal_id'] = reduced['internal_id'].astype(int)
    reduced["midline_msk"] = reduced.apply(lambda x: transform_midlines(x['midline'], x['align_seq'], x['align_pdb']), axis=1)
    return reduced


def pdb_first_gap(seq, start):
    ''' Defines where the first gap was placed. '''
    pos = np.array([i for i, letter in enumerate(seq) if letter == '-'])
    idgap = np.where(np.array(pos)==start)[0]
    if len(idgap)>0:
        idgap = idgap[0] 
    else:
        idgap = -1
    return idgap, pos


def pdb_size_start_gap(seq, start):
    ''' Determines the exact size of the gap region. '''
    len_gap = 0
    idgap, pos = pdb_first_gap(seq, start)
    if idgap>-1:
        pos = pos[idgap:]
        lastgap = np.nonzero(np.diff(pos)!=1)
        if len(lastgap[0])>0:
            len_gap = lastgap[0][0]
        else:
            len_gap = len(pos)
    return len_gap


def pos_is_gap(seq, start):
    ''' Counts how many extra positions should be added in case the place where
    the IDR/Poly should start is a gap. '''
    idgap, _ = pdb_first_gap(seq, start)
    if idgap>-1:
        len_gap = pdb_size_start_gap(seq, start)
        start=start+len_gap+idgap
    return start


def capture_idr_gaps(blast_over, prefix):
    ''' Modifies the coordinates and related fields to consider the informations
    extracted by the counting functions executed previously.
    cnt = increment every step of the loop
    cntb = increment only if there's no AA before a gap (compensates the gaps before the IDR start not counted)
    cnta = increment when there's no gap to make sure it stops counting at the IDR size
    '''
    if (prefix=="idr"):
        linker_id = 'internal_id'
    else:
        linker_id = prefix+'_name'
    blast_regs = blast_over.loc[:, [linker_id, prefix+'_seq_align_start', prefix+'_aa', 
                                    'last_part_seq', 'over_sz', prefix+'_seq_diff']]
    blast_regs_t = blast_regs.set_index(linker_id).T
    blast_regs_t = blast_regs_t.to_dict(orient='list')
    for k,v in blast_regs_t.items():
        cnt, cntb, cnta, push = 0,0,0,0
        st=''
        # Dealing with PDB alignments in the middle of the IDR
        if v[4]<0:
            push = int(abs(v[4]))
            aa = v[1][push:]
        else:
            aa = v[1][:int(v[3])]
        for p in v[2]:
            cnt+=1
            if p!='-':
                cnta+=1
                st+=p
            else:
                if cnta==0:
                    cntb+=1
            aa_pos = st.find(aa,0)
            if (aa_pos!=-1) and (cnta>=int(v[3])):
                break
        if (aa_pos==-1):
            aa_pos=0
        # Had to add this check because the gaps may happen after some AAs and
        # whithout them we can't find the start of the AA (tostines scenario)
        if aa_pos!=0:
            if cntb==0:
                before_aa = v[2][:aa_pos+1]
                cntb = before_aa.count('-')
        blast_regs_t[k][0] = int(v[0]+cntb+aa_pos)
        blast_regs_t[k].append(int(cnt-cntb-aa_pos-1))
    cols=list(blast_regs.columns)+['over_sz_seq']
    blast_regs = pd.DataFrame.from_dict(blast_regs_t, orient='index').reset_index()
    blast_regs.set_axis(cols, axis=1, inplace=True)
    blast_regs = blast_regs.loc[:, [linker_id, prefix+'_seq_align_start', 'over_sz_seq']]
    blast_over = pd.merge(blast_over, blast_regs, how='left', on=linker_id, suffixes=("_x", ""))
    return blast_over


def extract_gap_seq(pdb_data, prefix, linker_id):
    ''' Gets the coordinates of the sequence side of the alignment. THe numeric
    positions are required to all extractions of the other sides of the alignments.
    We needed to account for gaps inserted before, in the middle and inside the 
    IDR/Poly regions.'''
    blast_over = pdb_data.loc[~pdb_data['internal_id'].isna(),
                              ['internal_id', 'align_seq', prefix+'_start', 'seq_start',
                               'over_sz', prefix+'_name', prefix+'_aa']]
    #blast_over2 = pdb_data.loc[pdb_data[prefix+'_name'].isin(['A0A024QZ33_1', 'A0A024R4E5_1', 'A0A075B6E2_1', 'A0A024RB15_1', 'A0A1B0GTK5_1','A0A087WSY3_3', 'Q14162_2', 'Q96GM8_2', 'A0A087WSX5_5']),['internal_id', prefix+'_start', 'seq_start', 'pdb_start', 'over_sz', prefix+'_name', prefix+'_aa', 'over_dir']]; blast_over=blast_over2
    blast_over[prefix+'_seq_diff'] = blast_over[prefix+'_start']-blast_over['seq_start']
    blast_over[prefix+'_seq_new_1_start'] = blast_over.apply(lambda x: x[prefix+'_seq_diff'] if (x[prefix+'_seq_diff']>0) else 0, axis=1)
    blast_over = blast_over.sort_values(by=[prefix+'_name', 'over_sz'], ascending=[True, False]).reset_index(drop=True)
    # This deals with the cases where the position where the IDR should be is a GAP
    # Tricky thing: It is required to extract the sequence's first part, not the second
    blast_over[prefix+'_seq_new_2_start'] = blast_over.apply(lambda x: pos_is_gap(x['align_seq'], x[prefix+'_seq_new_1_start']), axis=1)
    blast_over['class_gaps_seq_start'] = blast_over[prefix+'_seq_new_1_start']-blast_over[prefix+'_seq_new_2_start']!=0
    blast_over['first_part_seq'] = blast_over.apply(lambda x: x['align_seq'][:int(x[prefix+'_seq_new_2_start'])], axis=1)
    blast_over[prefix+'_prev_gaps_seq'] = blast_over.first_part_seq.str.count('-')
    blast_over[prefix+'_seq_align_start'] = blast_over[prefix+'_seq_new_1_start']+blast_over[prefix+'_prev_gaps_seq']
    blast_over['class_gaps_seq_before'] = blast_over[prefix+'_seq_new_2_start']-blast_over[prefix+'_seq_align_start']!=0
    blast_over['last_part_seq'] = blast_over.apply(lambda x: x['align_seq'][int(x[prefix+'_seq_align_start']):], axis=1)
    # Here we count the gaps inside to resize the aligned area adjust the IDR 
    # alignment start, removing extra gaps remaining from the first cut.
    blast_over = capture_idr_gaps(blast_over, prefix)
    blast_over['add_prev_gaps_seq'] = blast_over[prefix+'_seq_align_start_x']-blast_over[prefix+'_seq_align_start']
    blast_over['class_gaps_seq_start2'] = blast_over['add_prev_gaps_seq']!=0
    # I reloaded the first part to count the official gaps before the region
    blast_over['first_part_seq'] = blast_over.apply(lambda x: x['align_seq'][:int(x[prefix+'_seq_align_start'])], axis=1)
    blast_over[prefix+'_prev_gaps_seq'] = blast_over.first_part_seq.str.count('-')
    # Just reload the last part to check if any gap remains in the beginning of the aligned IDR
    blast_over['last_part_seq'] = blast_over.apply(lambda x: x['align_seq'][int(x[prefix+'_seq_align_start']):], axis=1)
    # Now we know exactly where the IDR aligment ends
    #blast_over[prefix+'_seq_align_start'] = blast_over[prefix+'_seq_align_start']-1
    blast_over[prefix+'_seq_align_end'] = blast_over[prefix+'_seq_align_start']+blast_over['over_sz_seq']+1
    blast_over['region_seq'] = blast_over.apply(lambda x: x['align_seq'][int(x[prefix+'_seq_align_start']):int(x[prefix+'_seq_align_end'])], axis=1)
    blast_over['region_gaps'] = blast_over.region_seq.str.count('-')
    # Capturing the sequence real start and end. This info is needed to check if 
    # the PolyXY is inside the original IDR-PDB alignment or not
    blast_over[prefix+'_seq_rel_start'] = blast_over.apply(lambda x: x[prefix+'_start'] if (x['seq_start']<=x[prefix+'_start']) else x['seq_start'], axis=1)

    blast_over_f = blast_over.loc[:, [linker_id, prefix+'_seq_align_start', 
                                      prefix+'_seq_align_end', prefix+'_prev_gaps_seq', 
                                      'region_seq', prefix+'_seq_diff', 'over_sz_seq', 
                                      'class_gaps_seq_start', 'class_gaps_seq_before', 
                                      'class_gaps_seq_start2', prefix+'_seq_rel_start']]
    return blast_over_f


def extract_2D_info(pdb_data, dssp_path, prefix):
    ''' Now that we have the coordinates from the sequence side, extracting the
    2D structure information from the DSSP data using the alignment info. '''
    # Loading 2D sequences
    with open(dssp_path, 'rb') as handle:
        dssp_dict = pickle.load(handle)
    # Getting just the area of interest from the 2D sequence
    pdb_data['ss_2Dpdb'] = pdb_data.apply(lambda x: dssp_dict[x['pdb_name']][0][int(x['pdb_start'])-1:int(x['pdb_end'])], axis=1)
    # Adding the GAPs generated in the PDB side of the alignment to the 2d Sequence
    pdb_data['ss_2Dpdb_gaps'] = pdb_data.apply(lambda x: insert_gaps_2D(x['align_pdb'], x['ss_2Dpdb'], 1), axis=1)
    pdb_data['ss_2Dall_gaps'] = pdb_data.apply(lambda x: insert_gaps_2D(x['align_seq'], x['ss_2Dpdb_gaps'], 2, "#", "#"), axis=1)
    pdb_data['ss_region_pdb_gaps'] = pdb_data.apply(lambda x: x['ss_2Dpdb_gaps'][int(x[prefix+'_seq_align_start']):int(x[prefix+'_seq_align_end'])], axis=1)
    pdb_data[prefix+'_pdb_gaps'] = pdb_data.ss_region_pdb_gaps.str.count('-')
    pdb_data['final_region_ss_comp'] = pdb_data.apply(lambda x: insert_gaps_2D(x['region_seq'], x['ss_region_pdb_gaps'], 2, "#", "#"), axis=1)
    pdb_data[prefix+'_region_seq_gaps'] = pdb_data.final_region_ss_comp.str.count('#')
    pdb_data['ss_region_gaps_sz'] = pdb_data.final_region_ss_comp.str.len()
    pdb_data[prefix+'_seq_rel_end'] = pdb_data[prefix+'_seq_rel_start']+pdb_data['over_sz_seq']-pdb_data[prefix+'_region_seq_gaps']
        
    # Creating the official align SEQ, PDB and SS alignments without SEQ gaps.
    # All functions called after this point need to use them
    pdb_data['final_align_pdb'] = pdb_data.apply(lambda x: insert_gaps_2D(x['align_seq'], x['align_pdb'], 2, "#", "#"), axis=1)
    pdb_data['final_align_pdb'] = pdb_data['final_align_pdb'].str.replace('#', '')
    pdb_data['final_align_seq'] = pdb_data.align_seq.str.replace('-', '')
    pdb_data['final_region_seq'] = pdb_data.region_seq.str.replace('-', '')
    pdb_data['final_align_ss'] = pdb_data.ss_2Dall_gaps.str.replace('#', '')
    pdb_data['ss_region_noSeq_gaps'] = pdb_data.final_region_ss_comp.str.replace('#', '')
    pdb_data[prefix+'_final_seq_align_start'] = pdb_data[prefix+'_seq_align_start']-pdb_data[prefix+'_prev_gaps_seq']
    pdb_data[prefix+'_final_seq_align_end'] = pdb_data[prefix+'_seq_align_end']-pdb_data[prefix+'_region_seq_gaps']-pdb_data[prefix+'_prev_gaps_seq']
    pdb_data['ss_final_over_sz'] = pdb_data.final_region_ss_comp.str.len()-pdb_data[prefix+'_region_seq_gaps']
    
    pdb_data = pdb_data.drop(['ss_2Dpdb', 'ss_2Dpdb_gaps'], axis=1)
    
    return pdb_data


def extract_ss_surround_child(pdb_data, dssp_path, prefix, linker_id, delim=0):
    ''' Extract the sequences and secondary structures of the Poly or other
    region inside the IDR considering the center of the region for analysis 
    purposes. Here we can have an extra option of surrounding the region by
    delimited size, which ignores the IDR region and used the fixed size.'''
    
    with open(dssp_path, 'rb') as handle:
        dssp_dict = pickle.load(handle)
    cols = ['internal_id', prefix+'_name', 'final_align_pdb', 'final_align_seq', 'final_align_ss', 
            prefix+'_final_seq_align_start', prefix+'_final_seq_align_end',
            'ss_final_over_sz', 'idr_final_seq_align_start', 'idr_final_seq_align_end']
    
    ss_over = pdb_data.loc[~pdb_data['internal_id'].isna(), cols]
    ss_over['ss_sz_left'] = (ss_over['ss_final_over_sz']*0.5).apply(np.floor)
    ss_over['ss_sz_right'] = (ss_over['ss_final_over_sz']*0.5).apply(np.ceil)
    ss_over['seq_idr_left'] = ss_over.apply(lambda x: x['final_align_seq'][int(x['idr_final_seq_align_start']):int(x[prefix+'_final_seq_align_start']+x['ss_sz_left'])], axis=1)
    ss_over['seq_idr_right'] = ss_over.apply(lambda x: x['final_align_seq'][int(x[prefix+'_final_seq_align_end']-x['ss_sz_right']):int(x['idr_final_seq_align_end'])], axis=1)
    ss_over['pdb_idr_left'] = ss_over.apply(lambda x: x['final_align_pdb'][int(x['idr_final_seq_align_start']):int(x[prefix+'_final_seq_align_start']+x['ss_sz_left'])], axis=1)
    ss_over['pdb_idr_right'] = ss_over.apply(lambda x: x['final_align_pdb'][int(x[prefix+'_final_seq_align_end']-x['ss_sz_right']):int(x['idr_final_seq_align_end'])], axis=1)
    ss_over['ss_idr_left'] = ss_over.apply(lambda x: x['final_align_ss'][int(x['idr_final_seq_align_start']):int(x[prefix+'_final_seq_align_start']+x['ss_sz_left'])], axis=1)
    ss_over['ss_idr_right'] = ss_over.apply(lambda x: x['final_align_ss'][int(x[prefix+'_final_seq_align_end']-x['ss_sz_right']):int(x['idr_final_seq_align_end'])], axis=1)
    
    # Here I limited to the area before the IDR limited by the IDR size.
    # This delimitation can give us some clue about some polyXY structure roles
    ss_over['seq_surr_idr_left'] = ss_over.apply(lambda x: x['seq_idr_left'][:-(int(x['ss_sz_left']))], axis=1)
    ss_over['seq_surr_idr_right'] = ss_over.apply(lambda x: x['seq_idr_right'][int(x['ss_sz_right']):], axis=1)
    ss_over['pdb_surr_idr_left'] = ss_over.apply(lambda x: x['pdb_idr_left'][:-(int(x['ss_sz_left']))], axis=1)
    ss_over['pdb_surr_idr_right'] = ss_over.apply(lambda x: x['pdb_idr_right'][int(x['ss_sz_right']):], axis=1)
    ss_over['ss_surr_idr_left'] = ss_over.apply(lambda x: x['ss_idr_left'][:-(int(x['ss_sz_left']))], axis=1)
    ss_over['ss_surr_idr_right'] = ss_over.apply(lambda x: x['ss_idr_right'][int(x['ss_sz_right']):], axis=1)
    
    if (delim!=0):
        
        ss_over['seq_align_left'] = ss_over.apply(lambda x: x['final_align_seq'][:int(x[prefix+'_final_seq_align_start']+x['ss_sz_left'])], axis=1)
        ss_over['seq_align_right'] = ss_over.apply(lambda x: x['final_align_seq'][int(x[prefix+'_final_seq_align_start']+x['ss_sz_left']):], axis=1)
        ss_over['pdb_align_left'] = ss_over.apply(lambda x: x['final_align_pdb'][:int(x[prefix+'_final_seq_align_start']+x['ss_sz_left'])], axis=1)
        ss_over['pdb_align_right'] = ss_over.apply(lambda x: x['final_align_pdb'][int(x[prefix+'_final_seq_align_start']+x['ss_sz_left']):], axis=1)
        ss_over['ss_align_left'] = ss_over.apply(lambda x: x['final_align_ss'][:int(x[prefix+'_final_seq_align_start']+x['ss_sz_left'])], axis=1)
        ss_over['ss_align_right'] = ss_over.apply(lambda x: x['final_align_ss'][int(x[prefix+'_final_seq_align_start']+x['ss_sz_left']):], axis=1)
        ss_over['diff_idr_left'] = delim-ss_over.pdb_align_left.str.len()
        ss_over['diff_idr_right'] = delim-ss_over.pdb_align_right.str.len()
        ss_over['seq_delim_left'] = ss_over.apply(lambda x: '|'*x['diff_idr_left']+x['seq_align_left'] if x['diff_idr_left']>0 else x['seq_align_left'][abs(x['diff_idr_left']):], axis=1)
        ss_over['seq_delim_right'] = ss_over.apply(lambda x: x['seq_align_right']+'|'*x['diff_idr_right'] if x['diff_idr_right']>0 else x['seq_align_right'][:delim], axis=1)
        ss_over['pdb_delim_left'] = ss_over.apply(lambda x: '|'*x['diff_idr_left']+x['pdb_align_left'] if x['diff_idr_left']>0 else x['pdb_align_left'][abs(x['diff_idr_left']):], axis=1)
        ss_over['pdb_delim_right'] = ss_over.apply(lambda x: x['pdb_align_right']+'|'*x['diff_idr_right'] if x['diff_idr_right']>0 else x['pdb_align_right'][:delim], axis=1)
        ss_over['ss_delim_left'] = ss_over.apply(lambda x: '|'*x['diff_idr_left']+x['ss_align_left'] if x['diff_idr_left']>0 else x['ss_align_left'][abs(x['diff_idr_left']):], axis=1)
        ss_over['ss_delim_right'] = ss_over.apply(lambda x: x['ss_align_right']+'|'*x['diff_idr_right'] if x['diff_idr_right']>0 else x['ss_align_right'][:delim], axis=1)
    
    ss_over_f = ss_over.loc[:, [linker_id, 'seq_idr_left', 'seq_idr_right',
                                'pdb_idr_left', 'pdb_idr_right', 'ss_idr_left', 
                                'ss_idr_right', 'seq_surr_idr_left',
                                'seq_surr_idr_right', 'pdb_surr_idr_left',
                                'pdb_surr_idr_right', 'ss_surr_idr_left', 
                                'ss_surr_idr_right', 'seq_delim_left', 
                                'seq_delim_right', 'pdb_delim_left', 
                                'pdb_delim_right', 'ss_delim_left', 'ss_delim_right']]
    return ss_over_f


def extract_ss_surround(pdb_data, dssp_path, prefix, linker_id):
    ''' Extract the sequences and secondary structures of the IDR considering 
    the center of the region for analysis purposes. '''
    # Loading 2D sequences
    with open(dssp_path, 'rb') as handle:
        dssp_dict = pickle.load(handle)
    
    cols = ['internal_id', prefix+'_name', 'final_align_pdb', 'final_align_seq', 'final_align_ss', 
            prefix+'_final_seq_align_start', prefix+'_final_seq_align_end', 'ss_final_over_sz']
    ss_over = pdb_data.loc[~pdb_data['internal_id'].isna(), cols]

    # Here we get the complete alignment region including the IDR
    ss_over['ss_sz_left'] = round(ss_over['ss_final_over_sz']*0.5, 0)
    ss_over['ss_sz_right'] = (ss_over['ss_final_over_sz']*0.5).apply(np.ceil)
    ss_over['seq_surr_left'] = ss_over.apply(lambda x: x['final_align_seq'][:int(x[prefix+'_final_seq_align_start']+x['ss_sz_left'])], axis=1)
    ss_over['seq_surr_right'] = ss_over.apply(lambda x: x['final_align_seq'][int(x[prefix+'_final_seq_align_end']-x['ss_sz_right']):], axis=1)
    ss_over['pdb_surr_left'] = ss_over.apply(lambda x: x['final_align_pdb'][:int(x[prefix+'_final_seq_align_start']+x['ss_sz_left'])], axis=1)
    ss_over['pdb_surr_right'] = ss_over.apply(lambda x: x['final_align_pdb'][int(x[prefix+'_final_seq_align_end']-x['ss_sz_right']):], axis=1)
    ss_over['ss_surr_left'] = ss_over.apply(lambda x: x['final_align_ss'][:int(x[prefix+'_final_seq_align_start']+x['ss_sz_left'])], axis=1)
    ss_over['ss_surr_right'] = ss_over.apply(lambda x: x['final_align_ss'][int(x[prefix+'_final_seq_align_end']-x['ss_sz_right']):], axis=1)
    
    ss_over_f = ss_over.loc[:, [linker_id, 
                                'seq_surr_left', 'seq_surr_right',
                                'pdb_surr_left', 'pdb_surr_right',
                                'ss_surr_left', 'ss_surr_right']]
    return ss_over_f


def find_first_aa(seq):
    ''' Get the position of the first gap in the PDB side. '''
    index_aa = re.search("\w+", seq)
    if index_aa:
        start_gap = index_aa.start()
    else:
        start_gap = len(seq)
    return int(start_gap)


def get_align_coords(part_seq, full_seq, internal_id):
    ''' Get the final start and end coordinates of the aligned region. '''
    align_start = full_seq.find(part_seq)
    align_end = align_start+len(part_seq)-1
    return pd.Series((v for v in [align_start, align_end]))


def extract_gap_pdb(pdb_data, prefix, linker_id):
    ''' Gets the coordinates of the PDB side of the alignment and final coordinates
    required to make sure the correct AUTH value is extracted to adjust the final 
    region coordinates. '''
    blast_over = pdb_data.loc[~pdb_data['internal_id'].isna(),
                              ['internal_id', prefix+'_seq_align_start', 'align_pdb', 
                               'pdb_name', 'pdb_start', 'over_sz', 'over_sz_seq', 
                               prefix+'_name', prefix+'_pdb_gaps', 'ss_region_gaps_sz', 'final_region_ss_comp', 
                               'pdb_aa', 'pdb_end', prefix+'_seq_align_end']]
    
    # If the region starts with GAPs, counts it and subtract from the sequence start
    blast_over[prefix+'_pdb_new_end'] = blast_over[prefix+'_seq_align_start']+blast_over['ss_region_gaps_sz']
    blast_over['region_pdb'] = blast_over.apply(lambda x: x['align_pdb'][int(x[prefix+'_seq_align_start']):int(x[prefix+'_pdb_new_end'])], axis=1)
    blast_over['pdb_inicial_gaps'] = blast_over.apply(lambda x: find_first_aa(x['region_pdb']), axis=1)
    blast_over[prefix+'_pdb_new_2_start'] = blast_over[prefix+'_seq_align_start']#-blast_over['pdb_inicial_gaps']
    blast_over['class_gaps_pdb_start'] = blast_over[prefix+'_pdb_new_2_start']-blast_over[prefix+'_seq_align_start']!=0
    
    # I also need the gaps before the IDR in the PDB sequence
    blast_over['first_part_pdb'] = blast_over.apply(lambda x: x['align_pdb'][:int(x[prefix+'_pdb_new_2_start'])], axis=1)
    blast_over['prev_gaps_pdb'] = blast_over.first_part_pdb.str.count('-')
            
    # This is the start for the original sequence. I need it to get the AA region
    # of the PDB side to count GAPs, but you'll realize that I use the seq_new_3_start
    # to find out the relative start in PDB, accounting for PDB and SEQ gaps
    blast_over[prefix+'_pdb_rel_start'] = blast_over[prefix+'_pdb_new_2_start']+blast_over['pdb_start']-blast_over['prev_gaps_pdb']
    blast_over['class_gaps_pdb_before'] = blast_over['prev_gaps_pdb']!=0
    # Discount the gaps and make sure the end position of the original sequence
    # is the same as the start if all the IDR is composed of gaps.
    blast_over[prefix+'_pdb_end_gaps'] = blast_over[prefix+'_pdb_gaps']-blast_over['pdb_inicial_gaps']
    blast_over['class_gaps_inside'] = blast_over[prefix+'_pdb_end_gaps']!=0
    blast_over['ss_over_new_sz'] = blast_over['ss_region_gaps_sz']-blast_over['pdb_inicial_gaps']
    blast_over['class_gaps_only'] = blast_over['ss_over_new_sz']==0
    blast_over['ss_over_new_sz'] = blast_over['ss_over_new_sz']-blast_over[prefix+'_pdb_end_gaps']
    # Define if one position will be subtracted. I needed to do this because
    # when the IDR is just a GAP I want to make sure the start and end position
    # are the same, otherwhise by standard the end position is always -1
    blast_over['subtract_end'] = blast_over['ss_over_new_sz'].apply(lambda x: x if x==0 else 1)
    blast_over[prefix+'_pdb_rel_end'] = blast_over[prefix+'_pdb_rel_start']+blast_over['ss_over_new_sz']-blast_over['subtract_end']
    
    blast_over['final_region_pdb_comp'] = blast_over.apply(lambda x: insert_gaps_2D(x['final_region_ss_comp'], x['region_pdb'], 2, '#', '#'), axis=1)
    blast_over['final_region_pdb'] = blast_over.final_region_pdb_comp.str.replace('#', '')
    
    # This are validation steps. I will use the sequence coords to cut the PDB region 
    # and align it to the PDB sequence to get complete start and end coordinates 
    # to check the coordinates generated above
    # REPEATING FIRST STEPS MADE ON SS BUT WITH PDB SEQUENCE
    blast_over['valid_pdb_align'] =  blast_over.apply(lambda x: x['pdb_aa'][int(x['pdb_start']):int(x['pdb_end'])], axis=1)
    blast_over['valid_pdb_align_gaps'] = blast_over.apply(lambda x: insert_gaps_2D(x['align_pdb'], x['valid_pdb_align'], 1), axis=1)
    # At this point I don't add the sequence gaps because they overlap the PDB AAs (we don't want that)
    blast_over['valid_pdb_region_gaps'] = blast_over.apply(lambda x: x['valid_pdb_align_gaps'][int(x[prefix+'_seq_align_start']):int(x[prefix+'_seq_align_end'])], axis=1)
    # Now I cut the gaps to be able to search the region in the original sequence again
    blast_over['valid_pdb_region'] = blast_over.valid_pdb_region_gaps.str.replace('-', '')
    blast_over[[prefix+'_valid_pdb_start', prefix+'_valid_pdb_end']] = blast_over.apply(lambda x: get_align_coords(x['valid_pdb_region'], x['pdb_aa'], x['internal_id']), axis=1)
    
            
    blast_over_f = blast_over.loc[:, [linker_id, prefix+'_pdb_rel_start', 
                                      prefix+'_pdb_rel_end', 'final_region_pdb', 
                                      'region_pdb', 'final_region_pdb_comp',
                                      'class_gaps_pdb_start', 'class_gaps_pdb_before',
                                      'class_gaps_only', 'class_gaps_inside',]]
                                      #,'valid_pdb_region', prefix+'_valid_pdb_start', prefix+'_valid_pdb_end']]
    return blast_over_f


def get_AAcomposition(seq_aa):
    ''' Calculate the Physical-Chemical proportions of the IDR and the second
    most common group. '''
    idr_int = [AA_CODE_LIST.index(aa) for aa in seq_aa]
    idr_group_tots = [sum([el in groups for el in idr_int]) for groups in AA_GROUPS]
    return pd.Series((v for v in idr_group_tots))#, idr_group_props


def calc_hamming_distance(string1, string2): 
    ''' Calculate the hamming distance of the target region. '''
    distance = 0
    L = len(string1)
    for i in range(L):
        if string1[i] != string2[i]:
            distance += 1
    return distance/len(string1)


def calc_sum_of_pairs(seq1, seq2, matrix, gap_s, gap_e, gap = True):
    ''' Calculate the sum of pair of the target region. '''
    if len(seq1)>0:
        for A,B in zip(seq1, seq2):
            diag = ('-'==A) or ('-'==B)
            yield (gap_e if gap else gap_s) if diag else matrix[(A,B)]
            gap = diag


def count_dssp(dssp_pdbs, add_group=False):
    ''' Counts the occurence of each kind of structure inside the region to
    generate the proportion of coverage. '''
    struct=-1
    count_lst, prop_lst, map_lst = [], [], []
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


def prepare_ss_counts(df_ss, linker_id, dir_name="", cut="", prefix="idr"):
    ''' Prepare the required data to extract the counts the occurence of each 
    kind of structure. '''
    if (dir_name!=""):
        dir_name = "_"+dir_name
    if (cut!=""):
        cut = "_"+cut
    cols = list(df_ss.columns)
    idr_ss = dict(zip(df_ss[cols[0]], df_ss[cols[1]]))
    if (cut=='_delim'):
        ss_counts, ss_props, ss_map = count_dssp(idr_ss, True)
    else:
        ss_counts, ss_props, _  = count_dssp(idr_ss)
        ss_map=np.empty([1,1])
    
    df_ss["str1"] = df_ss[cols[3]].str.replace("|", "")
    df_ss["str2"] = df_ss[cols[2]].str.replace("|", "")
    df_ss["hamm_dist"] = df_ss.apply(lambda x: calc_hamming_distance(x["str1"], x["str2"]) if len(x["str1"])>0 else np.nan, axis=1)
    df_ss['sum_of_pairs'] = df_ss.apply(lambda x: sum(calc_sum_of_pairs(x["str1"], x["str2"], blosum, -5, -1)) if len(x["str1"])>0 else np.nan, axis=1)
    
    # Merging the counts to add to the main dataframe
    ids_idrs = df_ss[linker_id].values
    ss_region = df_ss.iloc[:, 1].values
    pdb_region = df_ss.iloc[:, 2].values
    seq_region = df_ss.iloc[:, 3].values
    hamm_vals = df_ss.iloc[:, 6].values
    sop_vals = df_ss.iloc[:, 7].values
    cols = [linker_id, prefix+'_seq'+cut+dir_name, prefix+'_pdb'+cut+dir_name, 
            prefix+'_ss'+cut+dir_name, prefix+'_hammDist'+cut+dir_name, prefix+'_sumOfPairs'+cut+dir_name,
            'cnt_helix'+cut+dir_name, 'cnt_sheet'+cut+dir_name, 'cnt_coil'+cut+dir_name, 
            'cnt_unmodeled'+cut+dir_name, 'cnt_unfolded'+cut+dir_name, 
            'cnt_gaps'+cut+dir_name, 'prop_helix'+cut+dir_name, 'prop_sheet'+cut+dir_name, 
            'prop_coil'+cut+dir_name, 'prop_unmodeled'+cut+dir_name, 
            'prop_unfolded'+cut+dir_name, 'prop_gaps'+cut+dir_name]
    if (cut=='_delim'):
        cols.insert(12, 'cnt_noStruct'+cut+dir_name)
        cols.append('prop_noStruct'+cut+dir_name)
    df_dist2D = pd.DataFrame(np.hstack((ids_idrs[:,None], seq_region[:,None], 
                                        pdb_region[:,None], ss_region[:,None],
                                        hamm_vals[:,None], sop_vals[:,None], 
                                        ss_counts, ss_props)), columns=cols)
    df_dist2D = df_dist2D.astype({'cnt_helix'+cut+dir_name: 'int32', 
                                  'cnt_sheet'+cut+dir_name: 'int32',
                                  'cnt_coil'+cut+dir_name: 'int32', 
                                  'cnt_unfolded'+cut+dir_name: 'int32',
                                  'cnt_unmodeled'+cut+dir_name: 'int32', 
                                  'cnt_gaps'+cut+dir_name: 'int32', 
                                  'prop_helix'+cut+dir_name: 'float32', 
                                  'prop_sheet'+cut+dir_name: 'float32',
                                  'prop_coil'+cut+dir_name: 'float32', 
                                  'prop_unfolded'+cut+dir_name: 'float32', 
                                  'prop_unmodeled'+cut+dir_name: 'float32', 
                                  'prop_gaps'+cut+dir_name: 'float32'})
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


def get_dssp_idr(dssp_path, blast_over_path, df_idr_details, prefix="idr", delim=0):
    ''' Take all the cases that paired to PDB structures and extract the region
    aligned to the IDR related to 2D sequences, considering the gaps generated
    by the alignment. '''
    # Getting just the fields needed to extract SEQ and PDB positions
    
    start_time = time.time()
    if (prefix!="idr"):
        linker_id = prefix+"_name"
    else:
        linker_id = "internal_id"
    cols = ['internal_id', 'seq_name','pdb_name',  prefix+'_name', 'over_sz', 
            prefix+'_start', 'seq_start', 'pdb_start', 'pdb_end', 'seq_aa', 
            'pdb_aa', prefix+'_aa', 'over_dir', 'pdb_size']
    if (prefix!="idr"):
        cols.append('idr_final_seq_align_start')
        cols.append('idr_final_seq_align_end')
    
    pdb_data = df_idr_details.loc[~df_idr_details['pdb_name'].isna(), cols].reset_index(drop=True)
    
    if (prefix=="idr"):
        pdb_data = pdb_data.sort_values(by=['internal_id', prefix+'_name'])
    else:
        pdb_data = pdb_data.sort_values(by=prefix+'_name')
    pdb_data_ids = pdb_data.index.tolist()
    pdb_data = pdb_data.reset_index(drop=True)
    # Getting the aligned sequences from the disk
    blast_over = get_blast_seqs(df_idr_details, blast_over_path)
    
    # The duplicates happen when I use the internalID from the IDR in another
    # set of data. Removing the duplicates is enough
    if (prefix!="idr"):
        blast_over = blast_over.drop_duplicates()
        
    print("PDB alignment data extracted. Time: {0} seconds".format((time.time() - start_time)/60))
        
    pdb_data = pd.merge(pdb_data, blast_over, how='left', on='internal_id')
    # Calulating and extracting the SEQ side of the alignment considering the GAPs
    blast_det = extract_gap_seq(pdb_data, prefix, linker_id)
    pdb_data = pd.merge(pdb_data, blast_det, how='left', on=linker_id)
    # Extract 2D sequence based in the coordinates of the SEQ with GAPs
    pdb_data = extract_2D_info(pdb_data, dssp_path, prefix)
    # JUST TO HELP IN THE VALIDATION PROCESS
    if (prefix=="idr"):
        ss_over = extract_ss_surround(pdb_data, dssp_path, prefix, linker_id)
    else:
        ss_over = extract_ss_surround_child(pdb_data, dssp_path, prefix, linker_id, delim)
    pdb_data = pd.merge(pdb_data, ss_over, how='left', on=linker_id)
    
    # Calculating and extracting the PDB side of the alignment considering the GAPs
    blast_det = extract_gap_pdb(pdb_data, prefix, linker_id)
    pdb_data = pd.merge(pdb_data, blast_det, how='left', on=linker_id)
    pdb_data['region_seq_nogaps'] = pdb_data.region_seq.str.replace('-', '')
    
    # Counting the polarity over the aligned SEQ region
    pdb_data[['tot_polar_pdb','tot_non_polar_pdb']] = pdb_data.apply(lambda x: get_AAcomposition(x['region_seq_nogaps']), axis=1)
    
    print("All data collected. Generating the SS couting. Time: {0} seconds".format((time.time() - start_time)/60))

    # Counting the 2D structures in the IDR region
    df_ss = pdb_data.loc[:, [linker_id, 'ss_region_noSeq_gaps', 'final_region_pdb', 'final_region_seq']]
    df_2D_idr, _ = prepare_ss_counts(df_ss, linker_id, "", "", prefix)
    
    # Adding the midline to the ss file (not cutted, complete midlines)
    pdb_data_sel = pdb_data.loc[:, [linker_id, "midline", "midline_msk"]]
    df_2D_idr = pd.merge(df_2D_idr, pdb_data_sel, how='left', on=linker_id)
    
    # Counting the 2D structures to the left of the IDR region
    cols=["ss_surr_left", "pdb_surr_left", "seq_surr_left"]
    if (prefix!="idr"):
        cols=["ss_idr_left", "pdb_idr_left", "seq_idr_left"]
    df_ss = pdb_data.loc[:, [linker_id]+cols]
    df_left_2D, _ = prepare_ss_counts(df_ss, linker_id, "left", "idr", prefix)
    df_2D_idr = pd.merge(df_2D_idr, df_left_2D, how='left', on=linker_id)
    
    # Counting the 2D structures to the right of the IDR region
    cols=["ss_surr_right", "pdb_surr_right", "seq_surr_right"]
    if (prefix!="idr"):
        cols=["ss_idr_right", "pdb_idr_right", "seq_idr_right"]
    df_ss = pdb_data.loc[:, [linker_id]+cols]
    df_right_2D, _ = prepare_ss_counts(df_ss, linker_id, "right", "idr", prefix)
    df_2D_idr = pd.merge(df_2D_idr, df_right_2D, how='left', on=linker_id)
    
    if (prefix!="idr"):
        # Counting the 2D structures to the left of the IDR region
        df_ss = pdb_data.loc[:, [linker_id, 'ss_surr_idr_left', 'pdb_surr_idr_left', 'seq_surr_idr_left']]
        df_left_2D, _ = prepare_ss_counts(df_ss, linker_id, "left", "surr", prefix)
        df_2D_idr = pd.merge(df_2D_idr, df_left_2D, how='left', on=linker_id)
        
        # Counting the 2D structures to the right of the IDR region
        df_ss = pdb_data.loc[:, [linker_id, 'ss_surr_idr_right', 'pdb_surr_idr_right', 'seq_surr_idr_right']]
        df_right_2D, _ = prepare_ss_counts(df_ss, linker_id, "right", "surr", prefix)
        df_2D_idr = pd.merge(df_2D_idr, df_right_2D, how='left', on=linker_id)
        
        # Counting the 2D structures for the delimited region
        if (delim!=0):
            df_ss = pdb_data.loc[:, [linker_id, 'ss_delim_left', 'pdb_delim_left', 'seq_delim_left']]
            df_left_2D, ss_map_left = prepare_ss_counts(df_ss, linker_id, "left", "delim", prefix)
            df_2D_idr = pd.merge(df_2D_idr, df_left_2D, how='left', on=linker_id)
            
            df_ss = pdb_data.loc[:, [linker_id, 'ss_delim_right', 'pdb_delim_right', 'seq_delim_right']]
            df_right_2D, ss_map_right = prepare_ss_counts(df_ss, linker_id, "right", "delim", prefix)
            df_2D_idr = pd.merge(df_2D_idr, df_right_2D, how='left', on=linker_id)
            
            map_50aa = np.hstack([ss_map_left, ss_map_right])
            map_50aa = map_50aa[pdb_data_ids, :]
            np.savetxt(BASIS_PATH+prefix+"_map_50aa.csv", map_50aa, fmt="%d", delimiter=",")
              
    sel_cols = ['internal_id', prefix+'_name', 'ss_final_over_sz', 'tot_polar_pdb','tot_non_polar_pdb',
                prefix+'_pdb_rel_start', prefix+'_pdb_rel_end', 
                prefix+'_seq_rel_start', prefix+'_seq_rel_end', 
                prefix+'_seq_align_start', prefix+'_seq_align_end',
                prefix+'_final_seq_align_start', prefix+'_final_seq_align_end']
                #, 'valid_pdb_region', prefix+'valid_pdb_start', prefix+'valid_pdb_end']
    if (prefix=="idr"):
        sel_cols = sel_cols+['final_align_ss']
    # merging the PDB and Seq original coodinates to the dataframe
    pdb_coords = pdb_data.loc[:, sel_cols]
    pdb_coords = pdb_coords.astype({prefix+'_pdb_rel_start': 'int32', 
                                  prefix+'_pdb_rel_end': 'int32',
                                  prefix+'_seq_rel_start': 'int32',
                                  prefix+'_seq_rel_end': 'int32',
                                  prefix+'_seq_align_start': 'int32',
                                  prefix+'_seq_align_end': 'int32',
                                  prefix+'_final_seq_align_start': 'int32',
                                  prefix+'_final_seq_align_end': 'int32'})
                                  #prefix+'_valid_pdb_start': 'int32',
                                  #prefix+'_valid_pdb_end': 'int32'})
        
    df_idr_details = pd.merge(df_idr_details, pdb_coords, how='left', on=['internal_id', prefix+'_name'])
    
    # Adding real percentage over structure considering the structure counts
    df_2D_idr["over_ss_real_sz"] = (df_2D_idr['cnt_helix']+df_2D_idr['cnt_sheet']+
                                    df_2D_idr['cnt_coil']+ df_2D_idr['cnt_unfolded'])
    df_2D_idr["over_ss_real_perc"] = df_2D_idr["over_ss_real_sz"]/df_2D_idr[prefix+"_seq"].str.len()
    
    df_real_ss = df_2D_idr.loc[:, [linker_id, "over_ss_real_sz", "over_ss_real_perc"]]
    df_idr_details = pd.merge(df_idr_details, df_real_ss, how='left', on=linker_id)
    df_2D_idr.drop(['over_ss_real_sz', 'over_ss_real_perc'], inplace=True, axis=1)
    
    print("All SS transformations done. Final time: {0} seconds".format((time.time() - start_time)/60))
    
    return df_idr_details, df_2D_idr, pdb_data

def remove_no_ss(df_2D_idr, cutoff):
    ''' Filter the regions composed by unmodeled residues according to the cutoff
    criteria. '''
    df_2D_idr_v2 = df_2D_idr.copy()
    df_2D_idr_v2['unmodeled'] = df_2D_idr_v2.idr_ss.str.len()-df_2D_idr_v2['cnt_unmodeled']
    df_2D_idr_v2['frac_unmodeled'] = df_2D_idr_v2['unmodeled']/df_2D_idr_v2.idr_ss.str.len()
    df_2D_idr_v2['keep1'] = (df_2D_idr_v2['unmodeled']>=cutoff[1])
    df_2D_idr_v2['keep2'] = (df_2D_idr_v2['frac_unmodeled']>=cutoff[0])
    df_2D_idr_v2['keep'] = df_2D_idr_v2['keep1']|df_2D_idr_v2['keep2']
    lst_ids = df_2D_idr_v2.loc[~df_2D_idr_v2['keep'], 'internal_id'].tolist()
    return lst_ids, df_2D_idr_v2


def final_filtering(df_idr_details, df_2D_idr, pdb_data, basis_path, prefix, cutoff=(.5,10)):
    ''' Remove all the cases composed by missing regions according to the 
    established cutoff criteria to then filter the best candidate based on the 
    selected ordering criteria. '''

    print("\nStarting the final filtering of best candidates...")
    
    df_idr_details_v2 = df_idr_details.copy()
    ids1 = df_idr_details_v2.loc[~df_idr_details_v2["pdb_name"].isna(), "idr_name"]
    print("\nBefore filtering unmodeled by the strucuture:")
    print('Rows PDB: {0} - Unique {1} IDRs'.format(len(ids1), len(np.unique(ids1))))
    
    # Removing alignments covered by unmodeled residues and gaps
    lst_ids, df_2D_idr_v2 = remove_no_ss(df_2D_idr, cutoff)
    df_idr_details_v2.loc[df_idr_details_v2['internal_id'].isin(lst_ids), "pdb_name":"ss_seq"] = np.nan
    df_idr_details_v2.loc[df_idr_details_v2['internal_id'].isin(lst_ids), "ss_final_over_sz":] = np.nan
   
    ids1 = df_idr_details_v2.loc[~df_idr_details_v2["pdb_name"].isna(), "idr_name"]
    print("\nAfter filtering unmodeled by the strucuture:")
    print('Rows PDB: {0} - Unique {1} IDRs'.format(len(ids1), len(np.unique(ids1))))
    
    df_idr_details_v2.to_csv(basis_path+'EVAL_data_int_all.csv', index=False)

    # Keeping just the highest e-value now that we removed the unmodeled
    df_idr_details_v2 = df_idr_details_v2.sort_values(by=[prefix+'_name', 'over_perc', 'pdb_evalue', 'pdb_resolution', 'bitscore'], ascending=[True, False, True, True, False]).drop_duplicates(subset=[prefix+'_name'])
    #df_idr_details_v2 = df_idr_details_v2.sort_values(by=[prefix+'_name', 'over_perc', 'pdb_evalue', 'bitscore'], ascending=[True, False, True, False]).drop_duplicates(subset=[prefix+'_name'])
    ids1 = df_idr_details_v2.loc[~df_idr_details_v2["pdb_name"].isna(), "idr_name"]
    print("\nKeeping just the top alignment e-value:")
    print('{0} PDBs from {1} unique IDRs of a total of {2} IDRs.'.format(len(ids1), len(np.unique(ids1)), len(df_idr_details_v2)))

    # Organizing and keeping just the 2D data for the selected alignemnt
    lst_ids = df_idr_details_v2.loc[~df_idr_details_v2.pdb_name.isna(), "internal_id"].tolist()
    df_2D_idr = df_2D_idr.loc[df_2D_idr["internal_id"].isin(lst_ids), :]
    df_ids = df_idr_details_v2.loc[~df_idr_details_v2.pdb_name.isna(), ["internal_id", prefix+"_name"]]
    df_2D_idr = pd.merge(df_2D_idr, df_ids, how='left', on='internal_id')
    df_2D_idr.drop('internal_id', axis=1, inplace=True)
    df_2D_idr = df_2D_idr.sort_values(by=[prefix+'_name']).reset_index(drop=True)
    col = df_2D_idr.pop(prefix+'_name')
    df_2D_idr.insert(0, prefix+'_name', col)
    
    save_name_all = 'data_all_'+prefix+'.csv'
    df_idr_details_v2 = mark_no_homologs(basis_path, df_idr_details_v2, save_name_all)
    save_name = 'data_ss_'+prefix+'.csv'
    df_2D_idr.to_csv(basis_path+save_name, index=False)
    
    return df_idr_details_v2, df_2D_idr, save_name_all


def correct_pdb_coords(pdb_coords, df_idr_details, idr_all_path, prefix="idr"):
    ''' Adjust PDB coordinates using the AUTH information extracted from CIF files. '''
    
    # Extracting only the columns of interest from both files
    cols = ["pdb_name", "dbrefs_start", "dbrefs_auth_start"]
    pdb_coords_sel = pdb_coords.loc[pdb_coords["pdb_ref_db"]!="PDB", cols]
    pdb_coords_sel = pdb_coords_sel.astype({"dbrefs_start":"int32", "dbrefs_auth_start":"int32"})
    pdb_coords_sel["dbrefs_real_start"] = pdb_coords_sel["dbrefs_auth_start"]-pdb_coords_sel["dbrefs_start"]
    cols = [prefix+"_name", "pdb_name", prefix+"_pdb_rel_start", "over_ss_real_sz"]
    df_idr_details_sel = df_idr_details.loc[~df_idr_details["pdb_name"].isna(), cols]
    df_idr_details_sel = pd.merge(df_idr_details_sel, pdb_coords_sel, how="left", on="pdb_name")
    
    # Generating final coordinates
    df_idr_details_sel[prefix+'_dbrefs_auth_start'] = df_idr_details_sel[prefix+"_pdb_rel_start"]+df_idr_details_sel["dbrefs_real_start"]
    df_idr_details_sel[prefix+'_dbrefs_auth_end'] = df_idr_details_sel[prefix+"_pdb_rel_start"] + df_idr_details_sel["dbrefs_real_start"]+df_idr_details_sel['over_ss_real_sz']-1
    df_idr_details_sel = df_idr_details_sel.drop(['pdb_name', prefix+'_pdb_rel_start', 'over_ss_real_sz'], axis=1)
       
    df_idr_details = pd.merge(df_idr_details, df_idr_details_sel, how="left", on=prefix+"_name")
    df_idr_details.to_csv(idr_all_path, index=False)
    return df_idr_details


def main_merge_idrPdb(idrs_path, pdb_path, file_names, pdb_det_path, cutoff):
    ''' Parts 1 and 2 - Main function IDR-PDB, filtering just synthetic strucures. '''
    # This process can take several minutes to run.
    
    basis_path = resources.get_dir(pdb_path)+"/"+resources.get_filename(pdb_path)+"_"
    
    # Not sure if the separators can be adjusted in the blast file
    sep1 = "|"
    sep2 = "|"
    # Dumping numpy objects to disk to save memory
    pickle_sz=10e+05
    
    print("\nMaking the XML file provided by blastP easier to my use...")    
    # Starting with the extraction of each blast pair from the xml file
    blast_pickle = extract_blast_pairs(pdb_path, pickle_sz, sep1, sep2)
    df_idr = pd.read_csv(idrs_path)
    
    # Replicating the PDBs and IDRs to extract Positions, E-value and unique IDSs an all other relevant alignment data
    print("\nTaking the starts and ends of the matches from the blastP output.")
    
    idx_million = positions2array(df_idr, blast_pickle, (0,8,9), ['idr_start', 'idr_end'], pickle_sz, 'starts_ends')
    
    print("\nTaking the rest of the thingies I need to use from blastP output.")
    get_another_data(df_idr, blast_pickle, idx_million, (0,1,2,5,6,4,10,11,12,7,14,15,16,8,9), 
                                   ['idr_name', 'idr_size'], pickle_sz, (1,1), file_names)
    
    if len(file_names)==5:
        file_names.insert(0, 'starts_ends')
    
    # We started running for less significant e-values too, but data proved not relevant
    # So we dropped the other e-values: 0.01 (maybe significant) and 10 (not significant).
    # Check the function to know details about the other outputs. They were created for validation purposes
    df_idr_new, _, _ = generate_all_dfs(basis_path, df_idr, None, 10e-05, "idr", cutoff, False, file_names)
    # Removes synthetic PDB structures from the dataframe with candidates and adds the experiment resolution
    df_idr_new = apply_pdb_data(df_idr_new, pdb_det_path)
    # IMPORTANT: You can run the same process passing df_idr_new as 2nd parameter 
    # to add less significant alignments.
    # If they are not relevant, you may consider add the significance paramenter to local blastP
    # to reduce the number of alignments and improve performance.
    eval_file = 'EVAL_data_noSeqs_idr.csv'
    df_idr_new = mark_no_homologs(basis_path, df_idr_new, eval_file)
    return basis_path+eval_file
    

def main_ss_annotation(idrs_path, pdb_mask_path, ss_file_path, dssp_path, idr_fasta_path, file_names, cutoff, save_dup=True):
    ''' Parts 3 - Main function IDR-PDB, getting the 2D structure and selecting 
    best candidates. '''
    
    basis_path = resources.get_dir(idrs_path)+"/"+resources.get_filename(idrs_path)+"_"
    
    if len(file_names)==5:
        file_names.insert(0, 'starts_ends')
    
    # Even when I tried to fix the mixed types problem by setting the type as 
    # object for pdb_name and kept_on before saving the dataframe the warning 
    # was showed. It seams its just thrown in the loading, so I simply used the 
    # suggested low memory clause.
    df_idr_new = pd.read_csv(idrs_path, low_memory=False)
    df_idr_details = append_pdb_details(pdb_mask_path, df_idr_new)
    # Important, the ss fasta is unsordered. ALWAYS use it through key selection from the dictionary.
    df_idr_details = append_ss_details(ss_file_path, df_idr_details)
    # ========================== Add seqs info ===============================
    # Need to add sequence details first to validate the real sequence start/end
    # position of the alignment later
    df_idr_details = append_seq_details(idr_fasta_path, df_idr_details)
    df_idr_details, df_2D_idr, pdb_data = get_dssp_idr(dssp_path, basis_path+file_names[5]+".pickle", df_idr_details)
    
    if save_dup:
        file_dup_idr = basis_path+'EVAL_data_dup_all.csv'
        df_idr_details.to_csv(file_dup_idr, index=False)
        # df_idr_details = pd.read_csv(file_dup_idr)
        file_dup_ss = basis_path+'EVAL_data_dup_ss.csv'
        df_2D_idr.to_csv(file_dup_ss, index=False)
        # df_2D_idr = pd.read_csv(file_dup_ss)
        file_dup_pdb = basis_path+'EVAL_data_dup_pdb.csv'
        pdb_data.to_csv(file_dup_pdb, index=False)
        # pdb_data = pd.read_csv(file_dup_pdb)
    
    df_idr_details, df_2D_idr, idr_all_path = final_filtering(df_idr_details, df_2D_idr, pdb_data, basis_path, "idr", cutoff)
    return basis_path+idr_all_path
    
    
def main_pos_files(idr_all_path, pdb_files, path_pdb_files):
    ''' Parts 4 - Main function IDR-PDB, getting the real PDB coordinates from
    the PDB files and adding to the detailed file. '''
    
    basis_path = resources.get_dir(idr_all_path)+"/"+resources.get_filename(idr_all_path)+"_"
    
    df_idr_details = pd.read_csv(idr_all_path, low_memory=False)
    cif_ids = np.unique(df_idr_details.loc[~df_idr_details["pdb_name"].isna(), "pdb_id"]).tolist()
    missing = pdbDssp.manage_missing(path_pdb_files, cif_ids, "-c", False) # ?? Need to delete previous files ??
    if len(missing)>0:
        file_name = basis_path +"sitll_missing_cif_list.txt"
        resources.save_sep_file(missing, file_name)
        print("ATENTION: There are still some CIF files we were not able to download.\nCheck the file sitll_missing_cif_list.txt and try to download them manually.\nWe will move on now considering the files available on disk!")
    df_pdb_coords, df_chain = pdbDssp.extract_from_cif_all(path_pdb_files, cif_ids)
    path_coords = basis_path+"coords_pdb.csv"
    path_chains = basis_path+"chains_pdb.csv"
    df_pdb_coords.to_csv(path_coords, index=False)
    df_chain.to_csv(path_chains, index=False)
    correct_pdb_coords(df_pdb_coords, df_idr_details, idr_all_path, "idr")
    
    
    