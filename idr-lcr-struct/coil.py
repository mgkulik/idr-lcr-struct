#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 25 12:18:26 2021
tuna
@author: magoncal
"""

#import idr_pdb
from itertools import product
from Bio import SeqIO
import os
import numpy as np
import pandas as pd
import time

def accum(accmap, a, func=None, size=None, fill_value=0, dtype=None):
    """
    An accumulation function similar to Matlab's `accumarray` function.

    Parameters
    ----------
    accmap : ndarray
        This is the "accumulation map".  It maps input (i.e. indices into
        `a`) to their destination in the output array.  The first `a.ndim`
        dimensions of `accmap` must be the same as `a.shape`.  That is,
        `accmap.shape[:a.ndim]` must equal `a.shape`.  For example, if `a`
        has shape (15,4), then `accmap.shape[:2]` must equal (15,4).  In this
        case `accmap[i,j]` gives the index into the output array where
        element (i,j) of `a` is to be accumulated.  If the output is, say,
        a 2D, then `accmap` must have shape (15,4,2).  The value in the
        last dimension give indices into the output array. If the output is
        1D, then the shape of `accmap` can be either (15,4) or (15,4,1) 
    a : ndarray
        The input data to be accumulated.
    func : callable or None
        The accumulation function.  The function will be passed a list
        of values from `a` to be accumulated.
        If None, numpy.sum is assumed.
    size : ndarray or None
        The size of the output array.  If None, the size will be determined
        from `accmap`.
    fill_value : scalar
        The default value for elements of the output array. 
    dtype : numpy data type, or None
        The data type of the output array.  If None, the data type of
        `a` is used.

    Returns
    -------
    out : ndarray
        The accumulated results.

        The shape of `out` is `size` if `size` is given.  Otherwise the
        shape is determined by the (lexicographically) largest indices of
        the output found in `accmap`.


    Examples
    --------
    >>> from numpy import array, prod
    >>> a = array([[1,2,3],[4,-1,6],[-1,8,9]])
    >>> a
    array([[ 1,  2,  3],
           [ 4, -1,  6],
           [-1,  8,  9]])
    >>> # Sum the diagonals.
    >>> accmap = array([[0,1,2],[2,0,1],[1,2,0]])
    >>> s = accum(accmap, a)
    array([9, 7, 15])
    >>> # A 2D output, from sub-arrays with shapes and positions like this:
    >>> # [ (2,2) (2,1)]
    >>> # [ (1,2) (1,1)]
    >>> accmap = array([
            [[0,0],[0,0],[0,1]],
            [[0,0],[0,0],[0,1]],
            [[1,0],[1,0],[1,1]],
        ])
    >>> # Accumulate using a product.
    >>> accum(accmap, a, func=prod, dtype=float)
    array([[ -8.,  18.],
           [ -8.,   9.]])
    >>> # Same accmap, but create an array of lists of values.
    >>> accum(accmap, a, func=lambda x: x, dtype='O')
    array([[[1, 2, 4, -1], [3, 6]],
           [[-1, 8], [9]]], dtype=object)
    """

    # Check for bad arguments and handle the defaults.
    if accmap.shape[:a.ndim] != a.shape:
        raise ValueError("The initial dimensions of accmap must be the same as a.shape")
    if func is None:
        func = np.sum
    if dtype is None:
        dtype = a.dtype
    if accmap.shape == a.shape:
        accmap = np.expand_dims(accmap, -1)
    adims = tuple(range(a.ndim))
    if size is None:
        size = 1 + np.squeeze(np.apply_over_axes(np.max, accmap, axes=adims))
    size = np.atleast_1d(size)

    # Create an array of python lists of values.
    vals = np.empty(size, dtype='O')
    for s in product(*[range(k) for k in size]):
        vals[s] = []
    for s in product(*[range(k) for k in a.shape]):
        indx = tuple(accmap[s])
        val = a[s]
        vals[indx].append(val)

    # Create the output array.
    out = np.empty(size, dtype=dtype)
    for s in product(*[range(k) for k in size]):
        if vals[s] == []:
            out[s] = fill_value
        else:
            out[s] = func(vals[s])

    return out


def load_seqs(path, s_type='seq'):
    seq_dict = dict()
    for seq_record in SeqIO.parse(path, 'fasta'):
        name_seq = seq_record.id.split("|")
        if (s_type=='seq'):
            name_seq = name_seq[1]
        elif (s_type=='pdb'):
            name_seq = name_seq[1]+name_seq[2]
        seq_dict[name_seq] = [str(seq_record.seq), seq_record.description]
    return seq_dict


# def get_curated_scores(taboth_path):
#     '''Extract the coiled coils scores (initially) and try to get the predicted
#     scores over the curated regions.'''
#     oth_dict = dict()
#     with open(taboth_path, 'r') as handle:
#         for line in enumerate(handle):
#             rowsplit = line[1].rstrip("\n").split("\t")
#             uniprod_id = rowsplit[0]
#             if (uniprod_id not in oth_dict):
#                 oth_dict[uniprod_id] = ["", [], []]
#             if (rowsplit[2]!="") and (len(oth_dict[uniprod_id][1])==0):
#                 oth_dict[uniprod_id][1] = [float(x) for x in rowsplit[2].split(",")]
#             if (rowsplit[3]!=""):
#                  oth_dict[uniprod_id][2] = rowsplit[3:]
#             if (oth_dict[uniprod_id][0]!='coiled_coil_uniprot'):
#                 oth_dict[uniprod_id][0] = rowsplit[1]
#     return oth_dict


def extract_high_scores(scores, tsh=.5, kmer=20):
    start_end = []
    score_cnt=0
    for c in enumerate(scores):
        if (c[1]>tsh):
            score_cnt+=1
        else:
            if score_cnt >= kmer:
                start_end.append([c[0]+1-score_cnt, c[0]])
            score_cnt=0
    if score_cnt >= kmer:
        start_end.append([c[0]+2-score_cnt, c[0]+1])
    return start_end


def extract_pos_coils(uniprot, reg, name):
    reg_num = len(reg)
    rep_unip = np.repeat(uniprot, reg_num).reshape(reg_num, 1)
    oth_ids = np.arange(1, reg_num+1).reshape(reg_num, 1)
    rep_type = np.repeat(name, reg_num).reshape(reg_num, 1)
    regs_arr = np.concatenate((rep_unip, rep_type), axis=1)
    regs_arr = np.concatenate((regs_arr, np.array(reg)), axis=1)
    regs_arr = np.concatenate((regs_arr, oth_ids), axis=1).tolist()
    return regs_arr
    

def get_mobidb_other(taboth_path, seq_dict, tsh=.5, kmer=20, last_id=0):
    oth_dict = dict()
    i=last_id+1
    with open(taboth_path, 'r') as handle:
        for line in enumerate(handle):
            seq_len=0
            rowsplit = line[1].rstrip("\n").split("\t")
            if (rowsplit[0] in seq_dict):
                seq_len = len(seq_dict[rowsplit[0]][0])
                if rowsplit[1]=="coiled_coil_uniprot":
                    name="coil_curat"
                else:
                    name="coil_fells"
                # Using the scores to get the region in the sequence
                if (rowsplit[2]!=""):
                    scores = [float(x) for x in rowsplit[2].split(",")]
                    scores = scores[0:seq_len]
                    reg = extract_high_scores(scores, tsh, kmer)
                # Getting the regions for the curated coils
                if rowsplit[3]!="":
                    reg = [[int(ri) for ri in r.split("-")] for r in rowsplit[3:]]
                if len(reg)>0:
                    regs_arr = extract_pos_coils(rowsplit[0], reg, name)
                    for j in range(0, len(regs_arr)):
                        oth_dict[i+j] = regs_arr[j]
                    i=i+len(regs_arr)
    return oth_dict, i-1

    
# def create_slide_scores(scores, k_mer=14):
#     n = len(scores)
#     xcol = np.arange(n-k_mer+1)
#     inds_col = np.tile(np.arange(k_mer), (n-k_mer+1,1))
#     inds_row = np.tile(xcol, (k_mer,1)).transpose()
#     scores_slide = scores[inds_col + inds_row]
#     mn_scores = np.max(scores_slide, axis=1)
    

# def generate_coil_threshold(oth_dict):
#     cur_coil = list()
#     for key, val in oth_dict.items():
#         if (len(val[2])>0) and (len(val[1])>0):
#             i=0
#             for reg in val[2]:
#                 i+=1
#                 pos = reg.split("-")
#                 coil_start=int(pos[0])
#                 coil_end=int(pos[1])
#                 coil_size=coil_end-coil_start
#                 reg_name = key+"_"+str(i)
#                 scores = np.array(val[1][coil_start-1:coil_end-1])
#                 mn_scores = np.mean(scores)
#                 md_scores = np.median(scores)
#                 #tss_scores = (scores-mn_scores)**2
#                 pr50_scores = sum(scores>.5)/coil_size
#                 gm_scores = gmean(scores)
#                 cur_coil.append([reg_name, coil_start, coil_end, coil_size, scores, mn_scores, md_scores, pr50_scores, gm_scores])
#     cols = ["coil_name", "coil_start", "coil_end", "coil_size", "coil_scores", "mean_coil", "med_coil", "prop50_coil", "gmean_coil"]
#     pd_cur_coil = pd.DataFrame(cur_coil, columns=cols)
#     pd_cur_coil.to_csv('analysis/coil_curated.csv', index=False)


def extract_coils_from_seq(coils_scores_path, coils_file, coils_name, last_id=0):
    oth_dict = dict()
    seqs_coils = load_seqs(coils_scores_path+coils_file)
    i=last_id+1
    for key, val in seqs_coils.items():
        seq_arr = np.array([v for v in val[0]])
        has_x = np.append(seq_arr=='x', False)
        has_x1 = np.insert(has_x[:-1], 0, False)
        reg_start = np.nonzero(has_x & ~has_x1)[0]+1
        reg_end = np.nonzero(~has_x & has_x1)[0]
        reg_num = len(reg_end)
        if reg_num>0:
            oth_ids = np.arange(1, reg_num+1).reshape(reg_num, 1)
            rep_unip = np.repeat(key, reg_num).reshape(reg_num, 1)
            rep_type = np.repeat(coils_name, reg_num).reshape(reg_num, 1)
            reg_pos = np.append(reg_start.reshape(reg_num,1), reg_end.reshape(reg_num,1), axis=1)
            regs_arr = np.append(rep_unip, rep_type, axis=1)
            regs_arr = np.append(regs_arr, reg_pos, axis=1)
            regs_arr = np.append(regs_arr, oth_ids, axis=1).tolist()
            for j in range(0, len(regs_arr)):
                oth_dict[i+j] = regs_arr[j]
            i=i+len(regs_arr)
    return oth_dict, i-1
    

# def append_seq_name(coils_in, seqs_lst, i):
#     sh = len(coils_in)
#     ids_rep = np.repeat(seqs_lst[i], sh).reshape(sh,1)
#     coils_in2 = np.append(ids_rep, np.array(coils_in).reshape(sh,1), axis=1)
#     return coils_in2


# def extract_coils_data(coils_scores_path, coils_file, seqs_lst):
#     with open(coils_scores_path+coils_file, 'r') as handle:
#         coils_lst, coils_in, coils_in3 = [],[],[]
#         i, j = 0, 0
#         for line in enumerate(handle):
#             rowsplit = line[1].rstrip("\n")
#             pos=int(rowsplit[0:4])
#             #aa=rowsplit[5]
#             #weigth=rowsplit[7]
#             #score=float(rowsplit[11:16])
#             perc=float(rowsplit[19:22])
#             if (pos==1) and (j>0):
#                 coils_in3.append(append_seq_name(coils_in, seqs_lst, i))
#                 coils_in = []
#                 i+=1
#             j+=1
#             coils_lst.append(perc)
#             coils_in.append(pos)
#             if i>50:
#                 break
#     coils_in3.append(append_seq_name(coils_in, seqs_lst, i))
#     coils_info = np.concatenate(coils_in3, axis=0)
#     coils_perc = np.array(coils_lst)
#     return coils_perc, coils_info


def merge_coil_dicts(dicts, f_name, path):
    new_dict = dicts[0].copy()
    for i in range(0, len(dicts)-1):
        if i!=0:
            new_dict = new_dict.copy()
        new_dict.update(dicts[i+1])
    resources.save_pickle(new_dict, f_name, path)
    return new_dict


def get_coil_ids(df_idr, oth_path, sz, path):
    ''' Generates a numpy array of unique IDs and other relevant data.'''
    un_ids = np.empty((sz,6), dtype=object)
    idr_sizes = np.empty(sz, dtype=int)
    i,k = 0,0
    IDR_old=''
    start_time = time.time()
    with open(oth_path, 'rb') as handle:
        while 1:
            try:
                dt_oth = pickle.load(handle)
                for key, val in dt_oth.items():
                    if (IDR_old) != val[0]:
                        filtered = df_idr[df_idr['seq_name']==val[0]]
                        ids_idr = filtered.loc[:, ['seq_name', 'idr_name']].values
                        idr_size = filtered.loc[:, 'idr_size'].values
                    
                    # Extracting and merging the IDs to create a unique identifier
                    internal_id = np.arange(k,k+len(filtered)).reshape(len(filtered),1)
                    oth_name = val[0]+"_"+str(val[4])+"_"+val[1]
                    oth_ids = np.tile([oth_name, val[1]], (len(filtered),1))
                    internal_id2 = np.repeat(i+1,len(filtered)).reshape(len(filtered),1)
                    ids = np.append(ids_idr, oth_ids, axis=1)
                    ids = np.append(ids, internal_id, axis=1)
                    ids = np.append(ids, internal_id2, axis=1)
                    un_ids[k:k+len(oth_ids),:] = ids
                    idr_sizes[k:k+len(oth_ids)] = idr_size
                    k=k+len(filtered)                           
                    IDR_old = val[0]
                    i+=1
                    if i%10e+03==0:
                        print(i)
            except EOFError:
                break
    print(i)
    resources.save_pickle(un_ids, 'unique_ids_coil.pickle', path)
    resources.save_pickle(idr_sizes, 'idr_sizes_coil.pickle', path)
    tot_in_sec = time.time() - start_time
    print("--- %s seconds ---" % (tot_in_sec))


def get_coils_group(coils_ids, coils_perc, df_idr, coil_type):
    coils_sel = coils_ids[:, 3]==coil_type
    coil_ids_sel = coils_ids[coils_sel, 1]
    coil_perc_sel = coils_perc[coils_sel]
    df_idr_ids = np.array(df_idr['idr_name'])
    _, _, int_sel = np.intersect1d(df_idr_ids, coil_ids_sel, return_indices=True, assume_unique=True)
    coil_perc_group = pd.DataFrame(
        np.column_stack((coil_ids_sel[int_sel], 
                         coil_perc_sel[int_sel])), 
        columns=["idr_name", coil_type])
    coil_perc_group = coil_perc_group.astype({coil_type: 'float'})
    df_idr = pd.merge(df_idr, coil_perc_group, how='left', on='idr_name')
    df_idr[coil_type] = df_idr[coil_type].fillna(0)
    return df_idr


def annotate_coils(df_idr, un_ids, over_bool, over_perc):
    coils_ids = un_ids[over_bool,:]
    coils_perc = over_perc[over_bool]
    coil_types = ['coil_coils28', 'coil_coils28w', 'coil_coils21', 'coil_coils21w']
    for names in coil_types:
        df_idr = get_coils_group(coils_ids, coils_perc, df_idr, names)
    return df_idr


def generate_df_coil(starts_ends, un_ids, orig_over_perc, over_perc, pos_over, over_bool):
    cols = ["seq_name", "idr_name", "coil_num", "coil_type", "id_merge1", "id_merge2", 
          "idr_start", "idr_end", "coil_start", "coil_end", "over_start", "over_end", 
          "coil_ind_perc", "coil_all_perc", "kept"]
    new_dtypes = {"id_merge1": int, "id_merge2": int, "idr_start": int, 
                  "idr_end": int, "coil_start": int, "coil_end": int, 
                  "over_start": int, "over_end": int, "coil_ind_perc": np.float64,
                  "coil_all_perc": np.float64, "kept": bool}
    df_coils_info = pd.DataFrame(np.column_stack((un_ids, starts_ends, 
                                          pos_over, orig_over_perc,
                                          over_perc, over_bool)), columns=cols)
    df_coils_info = df_coils_info.astype(new_dtypes)
    return df_coils_info


def merge_coil_coverage(un_ids, sz_over):
    _, un_idx, rep_idx = np.unique(un_ids[:,1]+un_ids[:,3], return_index=True, return_inverse=True)
    new_sz_over = np.bincount(un_idx[rep_idx], weights=sz_over)
    new_sz_over = new_sz_over[un_idx]
    return new_sz_over, un_idx, rep_idx


def get_coils_overlaps(df_idr, basis_path, new_path, cut_off=(.5,30)):
    starts_ends = cross_idr_pdb.open_pickle('starts_ends_coil.pickle', basis_path)
    un_ids = cross_idr_pdb.open_pickle('unique_ids_coil.pickle', basis_path)
    coil_sizes = cross_idr_pdb.open_pickle('idr_sizes_coil.pickle', basis_path)
    pos_over, sz_over = cross_idr_pdb.define_overlaps(starts_ends)
    # Calculates the coverage before the coils merge to use in the dataframe
    # with Coils information (validation purposes)
    orig_over_perc = cross_idr_pdb.get_overlap_size(sz_over, coil_sizes)
    # Now merge multiple coils over the same IDR
    sz_over, un_idx, rep_idx = merge_coil_coverage(un_ids, sz_over)
    coil_sizes = coil_sizes[un_idx]
    over_perc = cross_idr_pdb.get_overlap_size(sz_over, coil_sizes)
    over_bool = (over_perc>0)&((sz_over>cut_off[1])|(over_perc>cut_off[0]))
    # Gets just the first occurrence of each coil to get the correct coverage
    # Decidaded not to add details for coils (start-end) in the main IDR dataframe to save time
    # and the occurrences are not that relevant after all
    df_idr = annotate_coils(df_idr, un_ids[un_idx, :], over_bool, over_perc)
    df_idr.to_csv(new_path, index=False)
    # Replicating the boolean indexes to be able to extract all small portions 
    # of the coils in the same IDR (validation)
    over_bool_rep = over_bool[rep_idx]
    df_coils_info = generate_df_coil(starts_ends, un_ids, orig_over_perc, over_perc[rep_idx], pos_over, over_bool_rep)
    return df_idr, df_coils_info
    

def run_all(basis_path):
    coils28_w, last_id = extract_coils_from_seq(coils_scores_path, 'wf_coils28.out', 'coil_coils28w')
    nw_coils28, last_id = extract_coils_from_seq(coils_scores_path, 'nwf_coils28.out', 'coil_coils28', last_id)
    coils21_w, last_id = extract_coils_from_seq(coils_scores_path, 'wf_coils21.out', 'coil_coils21w')
    nw_coils21, last_id = extract_coils_from_seq(coils_scores_path, 'nwf_coils21.out', 'coil_coils21', last_id)
    coils_dict = merge_coil_dicts([coils28_w, nw_coils28, coils21_w, nw_coils21], 'coils_dict.pickle', basis_path)
