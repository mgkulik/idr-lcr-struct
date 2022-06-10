#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May 25 17:03:50 2020

Here I take the IDR regions and generate basic physical-chemical 
transformations to evaluate their characteristics. I also extract the 
complete sequences with IDRs to use in external tools and analyses.

05.06.2020 - Added the C-terminal 50 to split the data


@author: mgkulk
"""
from Bio import SeqIO
from operator import itemgetter
import numpy as np
import statistics as stats
import time
import pandas as pd
import os

import extract_coils_scores as coil_scores

pd.options.display.max_columns = 30

#'/home/magoncal/Documents/data/projects/idr_cook/'
#analysis_path = input("BE CAREFULL!! Set the complete path to the folder of the ANALYSIS files: ")
#basis_path = input("BE CAREFULL!! Set the complete path to the folder of the BASE files: ")

# tabidr_path = basis_path+'mobidb_UP000005640.tab'
# fastaname = basis_path+'uniprot-proteome_UP000005640.fasta'
# fastaout = basis_path+'mobidb_UP000005640_comp.fasta'
# idrs_path = analysis_path+'data_idrs.csv'
# error_path = basis_path+'mobidb_error.txt'

#tabidr_path = input("Path for the disordered positions (.tab file): ")
#fastaname = input("Fasta path for the complete proteom: ")
#fastaout =  input("Path to save the complete output fasta: ")
#idrs_path = input("Path to save the .csv file with IDR characteristics: ")
#error_path = input("Path to save errors from the IDR extraction: ")
#taboth_path = input("Path for the other info for IDRs (.tab file): ")

# Selenocysteine = 21, U, Sec, nonpolar
# Pyrrolysine = 22, O, Pyl, polar
AA_CODE_LIST = ['?','A','R','N','D','C','Q','E','G','H','I','L','K','M','F','P','S','T','W','Y','V','U','O','B','Z','X']
AA_CODE_ABRV = ['?', 'Ala', 'Arg', 'Asn', 'Asp', 'Cys', 'Gln', 'Glu', 'Gly', 'His', 'Ile', 'Leu', 'Lys', 'Met', 'Phe', 'Pro', 'Ser', 'Thr', 'Trp', 'Try', 'Val', 'Sec', 'Pyl', 'Asx', 'Glx', 'Xaa']
#AA_GROUPS = {'Non-Polar': [14,18,1,1lst_idrs5,20,11,10,13], 'Polar': [19,8,16,17,5,3,6], 'Basic': [12,2,9], 'Acidic': [4,7], 'Aromatics': [14,18,19]}
AA_GROUPS = [[19,8,16,17,5,3,6,22], [14,18,1,15,20,11,10,13,21], [12,2,9], [4,7]]
AA_GROUPS_NAMES = ['Polar Uncharged', 'Non-Polar', 'Polar Basic', 'Polar Acidic']

def get_seq_ints(seq):
    ''' Gets sequence AAs and generate their equivalent numeric sequence'''
    seq_int = np.array([AA_CODE_LIST.index(aa) for aa in seq])
    seq_by_aa = np.bincount(seq_int)
    seq_by_aa = np.pad(seq_by_aa, (0, len(AA_CODE_LIST) - len(seq_by_aa)), 'constant')
    seq_by_group = [sum(seq_by_aa[groups]) for groups in AA_GROUPS]
    return seq_by_group

def read_fasta(fastaname, seq_ints=False):
    ''' Reads the fasta file and store the sequences and their int equivalents. '''
    seqs_by_group = np.zeros((75777,4))
    fasta_data = dict()
    j=0
    seqs_len = 0
    for seq_record in SeqIO.parse(fastaname, 'fasta'):
        fasta_data[seq_record.id.split('|')[1]] = seq_record
        seq_len = len(str(seq_record.seq))
        if seq_ints:
            seqs_by_group[j, :] = get_seq_ints(str(seq_record.seq))
            seqs_len = seqs_len + len(str(seq_record.seq))
        j+=1
    return fasta_data, seqs_by_group, seqs_len
#del locals()['fasta_data']

                
def get_idr_AAcomposition(idr_aa):
    ''' Calculate the Physical-Chemical proportions of the IDR and the second
    most common group. '''
    idr_int = [AA_CODE_LIST.index(aa) for aa in idr_aa]
    idr_group_tots = [sum([el in groups for el in idr_int]) for groups in AA_GROUPS]
    idr_group_props = [i/len(idr_aa) for i in idr_group_tots]
    idr_tops_diff = sorted(idr_group_tots, reverse=True)
    idr_tops_diff = idr_tops_diff[0]-idr_tops_diff[1]
    idr_tops_diff_prop = sorted(idr_group_props, reverse=True)
    idr_tops_diff_prop = idr_tops_diff_prop[0]-idr_tops_diff_prop[1]
    idr_groups_var = np.var(idr_group_tots)
    idr_other_groups = sum(idr_group_tots[2:])
    idr_tops_idx = sorted(range(len(idr_group_tots)), key=lambda k: idr_group_tots[k], reverse=True)

    if (idr_tops_diff==0)&((idr_tops_idx[0]==0)&(idr_tops_idx[1]==1)):
        if (idr_other_groups>0):
            idr_cat1 = 'Polar Uncharged'
        else:
            idr_cat1 = 'Non-Polar'
    else:
        idr_cat1 = AA_GROUPS_NAMES[idr_group_props.index(max(idr_group_props))]
    idr_cat2 = AA_GROUPS_NAMES[idr_tops_idx[1]]
    return idr_group_tots, idr_group_props, [idr_cat1, idr_cat2, idr_tops_diff, idr_tops_diff_prop, idr_groups_var]


def get_idr_bins(idr_size, bin_size=30, bin_max=300):
    if idr_size>=bin_max:
        bin_start = bin_max
    elif idr_size==bin_size:
        bin_start = bin_size
    elif int(idr_size%bin_size) == 0:
        bin_start = int(bin_size*((idr_size/bin_size)-1))
    else:
        bin_start = bin_size*int(idr_size/bin_size)
    if bin_start == 0:
        bin_str = "(20-"+str(bin_start+(bin_size-1))+")"
    elif bin_start < 300:
        bin_str = "("+str(bin_start)+"-"+str(bin_start+(bin_size-1))+")"
    else:
        bin_str = str(bin_start)+"+"
    return bin_start, bin_str

#seq_idrs = rowsplit
#seq = str(seq_val.seq)
#seq_len = len(seq_val)
#idr = seq_idrs[start_pos:][0]
def get_idr_data(seq_idrs, seq, seq_len, start_pos, inter_size=50, perc_size=.7):
    ''' Loops over the IDR file with seq names and IDR positions separated
    by tab and gets general size and position information.'''
    all_idrs = list()
    error_seq = list()
    crit_size = inter_size*perc_size
    seq_name = seq_idrs[0]
    idr_scores = ""
    i=1
    if start_pos==3:
        scores = seq_idrs[2].split(',')
        idr_code = seq_idrs[1]
        if idr_code=='0':
            idr_type='predicted'
        elif idr_code=='1':
            idr_type='indicator'
        elif idr_code=='2':
            idr_type='curated'
    else:
        idr_type='predicted'
        scores=[]
    for idr in seq_idrs[start_pos:]:
        idrs = idr.split("-")
        idr_name = seq_name+'_'+str(i)
        idr_num = i
        idr_start = int(idrs[0])
        idr_rel_start = idr_start/seq_len
        idr_end = int(idrs[1])
        idr_rel_end = idr_end/seq_len
        idr_size = idr_end-idr_start+1
        idr_rel_size = idr_size/seq_len
        idr_bin, idr_bin_lbl = get_idr_bins(idr_size)
        if (idr_end<=seq_len):
            idr_aa = seq[idr_start-1:idr_end]
            if len(scores)>1:
                idr_scores = scores[idr_start-1:idr_end]
                num_scores = [float(s) for s in idr_scores]
                idr_mean_score = stats.mean(num_scores)
                idr_median_score = stats.median(num_scores)
                idr_scores = ", ".join(idr_scores)
            else:
                idr_mean_score, idr_median_score = 0, 0
                idr_scores = ''
            idr_group_tots, idr_group_props, idr_cats = get_idr_AAcomposition(idr_aa)
            #if (idr_size>=crit_size and (seq_len-idr_end)<=(inter_size-crit_size)):
            #    idr_last50=1
            lst_idrs = [seq_name, seq_len, idr_name, idr_num, idr_aa, idr_start, idr_rel_start, idr_end, idr_rel_end, idr_size, idr_rel_size, idr_bin, idr_bin_lbl, idr_type, idr_scores, idr_mean_score, idr_median_score] + idr_group_tots + idr_group_props + idr_cats
            all_idrs.append(lst_idrs)
            i+=1
        else:
            error_seq.append([idr_name+' - '+idr, 'IDR region outside the sequence.'])
    return all_idrs, error_seq


def save_fastas(seq_lst_comp, fastaout):
    ''' Saving filtered fasta files to use in other tasks. The first file is 
    complete, the others segregate de sequences in sequences which N AAs from
    the C-terminal belong to IDRs.'''    
    with open(fastaout, 'w') as handle:
      SeqIO.write(seq_lst_comp, handle, "fasta")

      
def save_file(lst, path):
    with open(path, 'a') as file:
        for data in lst:
          file.write(data)

          
def extract_idrs(tab_path, fasta_data, source, error_path=''):
    seq_lst, idrs_info, error_lst = [], [], []
    seq_idr_lens = 0
    if source=='mobidb_':
        start_pos = 3
    else:
        start_pos = 1
    with open(tab_path, 'r') as handle:
        i=0
        for line in enumerate(handle):
            rowsplit = line[1].rstrip("\n").split("\t")
            i+=1
            #if i==9:
            #    break
            if rowsplit[start_pos] != '':
                try:
                    seq_name = rowsplit[0].split("_")
                    seq_val = fasta_data[seq_name[0]]
                    # Getting IDR properties and adding to a DataFrame 
                    idr_details, error_seq = get_idr_data(rowsplit, str(seq_val.seq), len(seq_val), start_pos)
                    if len(idrs_info)==0:
                        idrs_info = idr_details
                    else:
                        idrs_info = idrs_info + idr_details
                    if len(error_seq)>0:
                        error_lst = error_lst + error_seq
                    # Getting sequence info to generate fastas from the sequences with IDRs
                    seq_lst.append(seq_val)
                    seq_idr_lens = seq_idr_lens+len(seq_val)
                except Exception as e:
                    error_lst.append([rowsplit[0], 'Sequence not available in the Uniprot file.'])
    if len(error_lst)>0:
        if error_path!="":
            with open(error_path, 'w') as outfile:
                for error in error_lst:
                    outfile.write(error[0]+' - '+str(error[1])+'\n')
        else:
            print("There are errors and the path to save them was not informed. Please check.")
    return seq_lst, seq_idr_lens, idrs_info, error_lst


def rename_idrs(pd_idrs):
    _, counts_seq = np.unique(pd_idrs['seq_name'], return_counts=True)
    rep_ids=[]
    for i in counts_seq:
        rep_ids=rep_ids+list(np.arange(1,i+1))
    pd_idrs['idr_num'] = rep_ids
    pd_idrs['idr_name'] = pd_idrs.seq_name.str.cat(pd_idrs['idr_num'].astype(str),sep="_")
    return pd_idrs


def generate_df_idrs(colnames, idrs_info, idrs_path, source, idr_min_sz=20):
    pd_idrs = pd.DataFrame(idrs_info, columns=colnames)
    pd_idrs = pd_idrs.sort_values(by=['seq_name', 'idr_start', 'idr_name'])
    pd_idrs = pd_idrs.loc[pd_idrs['idr_size']>=idr_min_sz, :]
    pd_idrs['idr_bin_lbl'] = pd_idrs["idr_bin_lbl"].apply(lambda x: '(300-'+str(max(pd_idrs['idr_size']))+')' if(str(x) == '300+') else x)
    if source=='mobidb_':
        pd_idrs = rename_idrs(pd_idrs)
    pd_idrs = pd_idrs.sort_values(by=['idr_name'])
    pd_idrs.to_csv(idrs_path, index=False)
    return pd_idrs


def remove_duplicate_seqs(seq_lst):
    new_seq_lst = list()
    _, idx_unique = list(np.unique([item.id for item in seq_lst], return_index=True))
    for i in idx_unique:
        new_seq_lst.append(seq_lst[i])
    return new_seq_lst


def set_threshold_str(tsh):
    return str(tsh).replace("0.", "").ljust(2,"0")


# Possible scores: {0.125, 0.25, 0.375, 0.5, 0.625, 0.75, 0.875, 1.0}
def get_score2regions(pd_idrs, tabidr_path, tsh=.5, kmer=20):
    tab_rows = list()
    sep = '\t'
    df_scores = pd_idrs.loc[pd_idrs["idr_scores"]!="", ["idr_name", "idr_scores", "idr_start"]]
    df_scores = df_scores.set_index('idr_name').T
    dict_scores = df_scores.to_dict(orient='list')
    for key, val in dict_scores.items():
        print(key)
        scores = [float(x) for x in val[0].split(",")]
        reg = coil_scores.extract_high_scores(scores, tsh, kmer)
        reg = [[r+val[1]-1 for r in s] for s in reg]
        regs = sep.join([str(re[0])+"-"+str(re[1]) for re in reg if re[1]-re[0] >= kmer])
        if regs!="":
            row_content = key+sep+regs+'\n'
            tab_rows.append(row_content)
    tabidr_path = tabidr_path.replace(".tab", "_"+set_threshold_str(tsh)+".tab")
    tabidr_path = tabidr_path.replace("mobidb", "mobidbnew")
    tabidr_path = tabidr_path.replace("dbase", "dbase/scores_sym")
    save_file(tab_rows, tabidr_path)
    return tabidr_path
       

def run_all():
    #pd_idrs = pd.read_csv('analysis/data_idrs.csv')
    # Getting the sequence and information about its relation to the last 50 AAs of
    # the sequence
    idr_min_sz=20
    fasta_data, seq_count_group, tot_lens = read_fasta(fastaname, False)
    source = os.path.basename(tabidr_path).split('_')[0]+"_"
    seq_lst, seq_idr_lens, idrs_info, error_lst = extract_idrs(tabidr_path, fasta_data, source, error_path)
    new_seq_lst = remove_duplicate_seqs(seq_lst)
    
    colnames = ['seq_name', 'seq_len', 'idr_name', 'idr_num', 'idr_aa', 'idr_start', 
                'idr_rel_start', 'idr_end', 'idr_rel_end', 'idr_size', 'idr_rel_size', 
                'idr_bin', 'idr_bin_lbl', 'idr_type', 'idr_scores', 'idr_mean_score', 
                'idr_median_score', 'tot_polar', 'tot_non_polar', 'tot_basic', 
                'tot_acidic', 'prop_polar', 'prop_non_polar', 'prop_basic', 
                'prop_acidic', 'idr_1st_cat', 'idr_2nd_cat', 'idr_tops_diff', 
                'idr_tops_diff_prop', 'idr_gp_variance']
    pd_idrs = generate_df_idrs(colnames, idrs_info, idrs_path, source, idr_min_sz)
    
    tsh=.6
    tabidr_path60 = get_score2regions(pd_idrs, tabidr_path, tsh)
    source = os.path.basename(tabidr_path60).split('_')[0]+"_"
    _, _, idrs_info60, _ = extract_idrs(tabidr_path60, fasta_data, source, max_len)
    idrs_path60 = idrs_path.replace(".csv", set_threshold_str(tsh)+".csv").replace("analysis", "analysis/scores_sym")
    pd_idrs60 = generate_df_idrs(colnames, idrs_info60, idrs_path60, source, idr_min_sz)
        
    # Sort descending the Max IDR size
    #ord_seq_lst = sorted(seq_lst, key = itemgetter(0), reverse=True)
    save_fastas(new_seq_lst, fastaout)
    #idx_idrs_max = np.argsort(-np.array(idrs_max))
        
    # Count AA occurrence in all sequences and in the ones with IDRs
    all_seqs_aa = np.sum(seq_count_group, axis=0, dtype=np.int32)
    allSeqs_group_prop = all_seqs_aa/tot_lens
    idr_seqs_aa = np.sum(seq_count_group[idx_count_group, :], axis=0, dtype=np.int32)
    idrSeqs_group_prop = idr_seqs_aa/seq_idr_lens