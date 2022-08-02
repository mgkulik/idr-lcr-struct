#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May 25 17:03:50 2020

Here I take the IDR regions and generate basic physical-chemical 
transformations to evaluate their characteristics. I also extract the 
complete sequences with IDRs to use in external tools and analyses.

@author: mgkulk
"""
from Bio import SeqIO
from operator import itemgetter
import numpy as np
import statistics as stats
import time
import pandas as pd
import json
import os

import resources
from localcider.sequenceParameters import SequenceParameters

pd.options.display.max_columns = 30

def get_entry_info(vals, name):
    ''' Returns basic info from mobidb predictions dictionary. '''
    pred_score = ""
    regions = vals['regions']
    if name=="mobidb":
        pred_score = vals['scores']
    return regions, pred_score

def extract_json(filename, key, name):
    ''' Reads the json proteome downloaded from mobidb, extract IDR predictions
    based on the consensus data and save a simple tab separated file to disk. '''
    
    sep = '\t'

    with open(filename, 'r') as handle:
        json_data = [json.loads(line) for line in handle]
    
    # Understanding the structure
    # The organism file presents only the consensus data ...
    dt_mobi = list()
    for a in range(len(json_data)):
        if (key in json_data[a]):
            seq_acc = json_data[a]['acc']
            seq_vals = json_data[a]['sequence']
            mobidb_vals = json_data[a][key]
            mobi_lst = ""
            reg_data, scores = get_entry_info(mobidb_vals, name)
            if name=="mobidb":
                scores = ','.join([str(i) for i in scores])
            mobi_lst = seq_acc+sep+scores 
            for k in range(len(reg_data)):
                mobi_lst += sep+str(reg_data[k][0])+"-"+str(reg_data[k][1])
            dt_mobi.append(mobi_lst+"\n")
    
    new_name = resources.gen_filename(filename, "mobidb", "idr", "tab")
    resources.save_file(dt_mobi, new_name)
    return (new_name)
    

def get_idr_bins(idr_size, bin_size=30, bin_max=300):
    ''' Set the bin size adjusting some label details. '''
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
#idr = seq_idrs[0]
def get_idr_data(seq_idrs, seq, seq_len, inter_size=50, perc_size=.7):
    ''' Loops over the IDR file with seq names and IDR positions separated
    by tab and gets general size and position information.'''
    all_idrs = list()
    error_seq = list()
    crit_size = inter_size*perc_size
    seq_name = seq_idrs[0]
    idr_scores = ""
    i=1
    scores=[]
    for idr in seq_idrs[2:]:
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
            idr_group_tots, idr_group_props, idr_cats = resources.get_idr_AAcomposition(idr_aa)
            #if (idr_size>=crit_size and (seq_len-idr_end)<=(inter_size-crit_size)):
            #    idr_last50=1
            lst_idrs = [seq_name, seq_len, idr_name, idr_num, idr_aa, idr_start, idr_rel_start, idr_end, idr_rel_end, idr_size, idr_rel_size, idr_bin, idr_bin_lbl, idr_scores, idr_mean_score, idr_median_score] + idr_group_tots + idr_group_props + idr_cats
            all_idrs.append(lst_idrs)
            i+=1
        else:
            error_seq.append([idr_name+' - '+idr, 'IDR region outside the sequence.'])
    return all_idrs, error_seq

          
def extract_idrs(tab_path, fasta_data, error_path=''):
    ''' Opens and get the IDR data and fasta data of the sequences of interest.
    Manages all I/O validations and save errors to disk. '''
    seq_lst, idrs_info, error_lst = [], [], []
    seq_idr_lens = 0
    with open(tab_path, 'r') as handle:
        i=0
        for line in enumerate(handle):
            rowsplit = line[1].rstrip("\n").split("\t")
            i+=1
            try:
                seq_name = rowsplit[0].split("_")
                seq_val = fasta_data[seq_name[0]]
                # Getting IDR properties and adding to a DataFrame 
                idr_details, error_seq = get_idr_data(rowsplit, str(seq_val.seq), len(seq_val))
                if len(idrs_info)==0:
                    idrs_info = idr_details
                else:
                    idrs_info = idrs_info + idr_details
                if len(error_seq)>0:
                    error_lst = error_lst + error_seq
                # Getting sequence info to generate fastas from the sequences with IDRs
                seq_lst.append(seq_val)
                seq_idr_lens = seq_idr_lens+len(seq_val[0])
            except Exception as e:
                # Errors happen when the fasta file does not have the IDs provided
                # by MobiDB. Two major reasons are obsolete entry or IDR annotated
                # in longer isoform. Decided to export evidence of the errors in 
                # disk and drop the IDRs.
                error_lst.append([rowsplit[0], 'Sequence not available in the Uniprot file.'])
    if len(error_lst)>0:
        if error_path!="":
            with open(error_path, 'w') as outfile:
                for error in error_lst:
                    outfile.write(error[0]+' - '+str(error[1])+'\n')
        else:
            print("There are errors and the path to save them was not informed. Please check.")
    return seq_lst, seq_idr_lens, idrs_info, error_lst


def generate_df_idrs(colnames, idrs_info, idrs_path, idr_min_sz=20):
    ''' Sorts and re-labels the bin max size using the biggest size of all IDRs.
    Saves the CSV file to disk.'''    
    pd_idrs = pd.DataFrame(idrs_info, columns=colnames)
    pd_idrs = pd_idrs.sort_values(by=['seq_name', 'idr_start', 'idr_name'])
    pd_idrs = pd_idrs.loc[pd_idrs['idr_size']>=idr_min_sz, :]
    pd_idrs['idr_bin_lbl'] = pd_idrs["idr_bin_lbl"].apply(lambda x: '(300-'+str(max(pd_idrs['idr_size']))+')' if(str(x) == '300+') else x)
    pd_idrs = pd_idrs.sort_values(by=['idr_name'])
    pd_idrs['proteome_id'] = resources.extract_proteome(idrs_path)
    pd_idrs.to_csv(idrs_path, index=False)
    return pd_idrs


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


def report_execution(start_time, name):
    ''' Return for how long this step runned. '''
    tot_in_sec = time.time() - start_time
    print("\nProp ", name, " calculated...")
    print("--- %s seconds ---" % (tot_in_sec))


def get_cider_props(idr_details_path, prefix):
    ''' Calculates IDR properties using cider. '''
    df_idr_details = pd.read_csv(idr_details_path, low_memory=False)
    start_time = time.time()
    print("\nStarting the properties calculation with cider.\nRelax, this can take some hours...")
    seqs_base = df_idr_details.loc[:, [prefix+'_name', prefix+'_aa']]
    seqs_base[prefix+'_WithXU'] = seqs_base[prefix+'_aa'].apply(lambda x: 0 if x.find('(?:X|U)')==-1 else x.find('(?:X|U)'))
    seqs_base[prefix+'_NoXU'] = seqs_base[prefix+'_aa'].str.replace('(?:X|U)', '', regex=True)
    seqs_base[prefix+'_par'] = seqs_base[prefix+'_NoXU'].apply(lambda x: SequenceParameters(x))
    seqs_base['fcr'] = seqs_base[prefix+'_par'].apply(lambda x: x.get_FCR())
    report_execution(start_time, list(seqs_base.columns)[-1])
    seqs_base['ncpr'] = seqs_base[prefix+'_par'].apply(lambda x: x.get_NCPR())
    report_execution(start_time, list(seqs_base.columns)[-1])
    seqs_base['isoPoint'] = seqs_base[prefix+'_par'].apply(lambda x: x.get_isoelectric_point())
    report_execution(start_time, list(seqs_base.columns)[-1])
    seqs_base['molWeight'] = seqs_base[prefix+'_par'].apply(lambda x: x.get_molecular_weight())
    report_execution(start_time, list(seqs_base.columns)[-1])
    seqs_base['fracNeg'] = seqs_base[prefix+'_par'].apply(lambda x: x.get_fraction_negative())
    report_execution(start_time, list(seqs_base.columns)[-1])
    seqs_base['fracPos'] = seqs_base[prefix+'_par'].apply(lambda x: x.get_fraction_positive())
    report_execution(start_time, list(seqs_base.columns)[-1])
    seqs_base['countNeg'] = seqs_base[prefix+'_par'].apply(lambda x: x.get_countNeg())
    report_execution(start_time, list(seqs_base.columns)[-1])
    seqs_base['countPos'] = seqs_base[prefix+'_par'].apply(lambda x: x.get_countPos())
    report_execution(start_time, list(seqs_base.columns)[-1])
    seqs_base['countNeut'] = seqs_base[prefix+'_par'].apply(lambda x: x.get_countNeut())
    report_execution(start_time, list(seqs_base.columns)[-1])
    seqs_base['promDisord'] = seqs_base[prefix+'_par'].apply(lambda x: x.get_fraction_disorder_promoting())
    report_execution(start_time, list(seqs_base.columns)[-1])
    seqs_base['kappa'] = seqs_base[prefix+'_par'].apply(lambda x: x.get_kappa())
    report_execution(start_time, list(seqs_base.columns)[-1])
    seqs_base['omega'] = seqs_base[prefix+'_par'].apply(lambda x: x.get_Omega())
    report_execution(start_time, list(seqs_base.columns)[-1])
    seqs_base['meanNCharge'] = seqs_base[prefix+'_par'].apply(lambda x: x.get_mean_net_charge())
    report_execution(start_time, list(seqs_base.columns)[-1])
    seqs_base['meanHydro'] = seqs_base[prefix+'_par'].apply(lambda x: x.get_mean_hydropathy())
    report_execution(start_time, list(seqs_base.columns)[-1])
    seqs_base['uverskyHydro'] = seqs_base[prefix+'_par'].apply(lambda x: x.get_uversky_hydropathy())
    report_execution(start_time, list(seqs_base.columns)[-1])
    seqs_base['ppiiPropHilser'] = seqs_base[prefix+'_par'].apply(lambda x: x.get_PPII_propensity(mode='hilser'))
    report_execution(start_time, list(seqs_base.columns)[-1])
    seqs_base['ppiiPropCreamer'] = seqs_base[prefix+'_par'].apply(lambda x: x.get_PPII_propensity(mode='creamer'))
    report_execution(start_time, list(seqs_base.columns)[-1])
    seqs_base['ppiiPropKallenbach'] = seqs_base[prefix+'_par'].apply(lambda x: x.get_PPII_propensity(mode='kallenbach'))
    report_execution(start_time, list(seqs_base.columns)[-1])
    seqs_base['delta'] = seqs_base[prefix+'_par'].apply(lambda x: x.get_delta())
    report_execution(start_time, list(seqs_base.columns)[-1])
    seqs_base['deltaMax'] = seqs_base[prefix+'_par'].apply(lambda x: x.get_deltaMax())
    report_execution(start_time, list(seqs_base.columns)[-1])
    
    seqpar_path = resources.get_dir(idr_details_path)+"/"+resources.get_filename(idr_details_path)+"_data_idr_properties.csv"
    seqs_base.to_csv(seqpar_path, index=False)
    tot_in_sec = time.time() - start_time
    print("\n--- Total cider time: %s seconds ---" % (tot_in_sec))


def run_all(fastaname, tabidr_path, use_toolscores=False):

    idr_min_sz=20
    fasta_data, seq_count_group, tot_lens = resources.read_fasta(fastaname)
    error_path = resources.gen_filename(tabidr_path, "mobidb", "idr", "err")
    seq_lst, seq_idr_lens, idrs_info, error_lst = extract_idrs(tabidr_path, fasta_data, error_path)
    
    if len(error_lst) > 0:
        print("There are errors, please check the file created on disk.")
    
    colnames = ['seq_name', 'seq_len', 'idr_name', 'idr_num', 'idr_aa', 'idr_start', 
                'idr_rel_start', 'idr_end', 'idr_rel_end', 'idr_size', 'idr_rel_size', 
                'idr_bin', 'idr_bin_lbl', 'idr_scores', 'idr_mean_score', 
                'idr_median_score', 'tot_polar', 'tot_non_polar', 'tot_basic', 
                'tot_acidic', 'prop_polar', 'prop_non_polar', 'prop_basic', 
                'prop_acidic', 'idr_1st_cat', 'idr_2nd_cat', 'idr_tops_diff', 
                'idr_tops_diff_prop', 'idr_gp_variance']
    idrs_path = resources.gen_filename(tabidr_path, "mobidb", "idr_details", "csv")
    pd_idrs = generate_df_idrs(colnames, idrs_info, idrs_path, idr_min_sz)
        
    if use_toolscores:
        # We end up not using this. The idea was to change the consensus percentage
        # and evaluate the curve compared to the confirmed IDRs.
        tsh=.6
        tabidr_path60 = get_score2regions(pd_idrs, tabidr_path, tsh)
        source = os.path.basename(tabidr_path60).split('_')[0]+"_"
        _, _, idrs_info60, _ = extract_idrs(tabidr_path60, fasta_data, source, max_len)
        idrs_path60 = idrs_path.replace(".csv", set_threshold_str(tsh)+".csv")
        pd_idrs60 = generate_df_idrs(colnames, idrs_info60, idrs_path60, idr_min_sz)
        

    fastaout = resources.gen_filename(tabidr_path, "mobidb", "idr", "fasta")
    resources.save_fastas(seq_lst, fastaout)
    #idx_idrs_max = np.argsort(-np.array(idrs_max))
     
    # Dropped: Compare the proportion of groups of AAs in the complete 
    # proteome and only in the sequences with IDRs. Kept commented because it 
    # uses some of the data generated as output of the read_fasta function.
    # Count AA occurrence in all sequences and in the ones with IDRs
    # all_seqs_aa = np.sum(seq_count_group, axis=0, dtype=np.int32)
    # allSeqs_group_prop = all_seqs_aa/tot_lens
    # idr_seqs_aa = np.sum(seq_count_group[idx_count_group, :], axis=0, dtype=np.int32)
    # idrSeqs_group_prop = idr_seqs_aa/seq_idr_lens
    
    return (idrs_path, fastaout)