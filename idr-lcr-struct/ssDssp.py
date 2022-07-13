#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 11 15:11:57 2022

@author: magoncal
"""

from Bio import SeqIO, pairwise2
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.PDB import PDBParser, MMCIFParser
from Bio.PDB.DSSP import DSSP
import pandas as pd
import numpy as np
import math
import re
import gzip, shutil, pickle, time

import os, subprocess, shlex
from datetime import datetime

import resources

path_seqres = "/home/magoncal/Documents/data/projects/poly_cook/2022_07_11_pdb_seqres.txt"
path_pdb_files = "/home/magoncal/Documents/data/projects/idr_cook/pdb_files_missing_ss"
path_masked = "/home/magoncal/Documents/data/projects/poly_cook/2021_12_06_pdb_seqres_masked.txt"


def download_seqs_file():
    "https://ftp.wwpdb.org/pub/pdb/derived_data/pdb_seqres.txt.gz"


def extract_seqres(path_seqres):
    ''' Extract the IDs and sequences from the seqRes fasta file provided by PDB 
    sequences download.
    Outputs:
        seqres_data: Dictionary with biopython SeqRecord data;
        seqres_keys : Numpy array of PDB IDs and chains using undescore as separator;
        seqres_desc: Dictionary with description of each PDB chain;
        seqres_fnb: Set of PDB IDs that have at least one molecule of type 
                    different from protein (Fold on bind cases).
    '''
    # Sample header: 104l_A mol:protein length:166  T4 LYSOZYME
    
    seqres_data = dict()
    seqres_keys = list()
    seqres_desc = dict()
    seqres_fnb = list()
    for seq_record in SeqIO.parse(path_seqres, "fasta"):
        seq_desc = seq_record.description.split(' ')
        seq_type = seq_desc[1].split(':')[1]
        seq_size = int(seq_desc[2].split(':')[1])
        comp = seq_record.id.split('_')
        if (seq_type=='protein'):
            comp = seq_record.id.split('_')
            name = comp[0].upper()+'_'+comp[1]
            seqres_data[name] = seq_record
            seqres_keys.append(name)
            seqres_desc[name] = ' '.join(seq_desc[4:]).lstrip()
        else:
            seqres_fnb.append(comp[0].upper())
    seqres_keys = np.array(seqres_keys)
    seqres_fnb = set(seqres_fnb)
    return seqres_data, seqres_keys, seqres_desc, seqres_fnb


def get_missing_names(seqres_keys, path_masked="", sep="_"):
    ''' Extract the list of files not available on the last extraction.
    If no path_masked is provided we assume this is the first execution, so all
    protein PDB files will be downloaded.
    '''
    if (path_masked!=""):
        ss_dict = resources.load_seqs(path_masked, "|", 'pdb')
        old_ss_keys = np.array(list(ss_dict.keys()))
        missing_names = np.setdiff1d(seqres_keys, old_ss_keys)
        missing_ids = np.unique([n.split(sep)[0] for n in missing_names]).tolist()
        obsolete_names = np.setdiff1d(old_ss_keys, seqres_keys)
    return missing_names, missing_ids, obsolete_names


def select_target(path_selected, un_pdb_ids, intersect=True, sep="_"):
    ''' Gets all files and extract just the target list based on an external 
    list of IDs. 
            un_pdb_ids must have the PDB ID WHITOUT CHAIN.
    '''
    files = sorted([f.path for f in os.scandir(path_selected) if (re.search(r'\W*.gz', f.path))])
    base_files = [os.path.splitext(os.path.splitext(os.path.basename(f))[0])[0] for f in files]
    sel_files = []
    if (intersect):
        for x, y in zip(base_files, files):
            if (x in un_pdb_ids):
                sel_files.append(y)
    else:
        sel_files = np.setdiff1d(un_pdb_ids, base_files)
    return sel_files


def manage_missing(path_pdb_files, missing_ids, file_name, del_prev=False):
    ''' Deletes old PDB/CIF files that showed errors previously so the script
    can try to run DSSP again and download the .pdb file using the bash file 
    provided by PDB and when not available the .cif file. Every 10 executions 
    it checks if it should abort in case some files remain missing to avoid 
    infinite loops. It will usually finish in the third execution when the 
    internet connection is OK.
    '''
    
    # Removing files already existing in the folder (There may be errors from 
    # last execution). It will make the execution slower. Avoid to do it daily.
    if del_prev:
        sel_files = select_target(path_pdb_files, missing_ids)
        if len(sel_files)>0 :
            try:
                for file in sel_files:
                    os.remove(file)
            except Exception as e:
                print([base, e])
        else:
            print("No files to delete.")
    
    # Now download all missing .pdb.gz files from PDB
    source_sh = os.path.join(path_pdb_files, "batch_download.sh")
    # Save missing files to disk, so the .sh provided by PDB can run over all 
    # downloads at once.
    missing = list(select_target(path_pdb_files, missing_ids, False))
    resources.save_sep_file(missing, file_name)
    prev_miss = len(missing)
    file_type = "-p"
    i = 0
    start_time = time.time()
    while len(missing) > 0:
        i += 1
        subprocess.call(shlex.split(f"{source_sh} -f {file_name} -o {path_pdb_files} {file_type}"))
        missing = list(select_target(path_pdb_files, missing, False))
        resources.save_sep_file(missing, file_name)
        if len(missing)!=prev_miss:
            prev_miss = len(missing)
        else:
            file_type = "-c"
        if (i%10==0):
            if (file_type=="-c" and len(missing)>0):
                print()
                print("{0} Attemps were already made and there are still missing PDB/CIF files.\nPlease download them manually.".format(str(i)))
                break
        print("{0} attempt(s) executed.".format(str(i)))
    print()
    tot_in_sec = time.time() - start_time
    print("--- %s seconds ---" % (tot_in_sec))
    os.remove(file_name)
    missing = list(select_target(path_pdb_files, missing_ids, False))
    return missing


def get_alignment(seq1, seq2, ma=2, nma=-1, og=-0.5, eg=-0.1):
    alignments = pairwise2.align.globalms(seq1, seq2, ma, nma, og, eg)
    if len(alignments[0].seqB) > len(seq1):
        alignments = pairwise2.align.globalms(seq1, seq2, ma, og, og, eg)
    return alignments[0].seqB


def get_ss_string(chain_part, base, seqres_data):
    new_lst = list()
    key = chain_part[0,1]
    chain_df = pd.DataFrame(chain_part[:,2:], index=chain_part[:,0], columns=['seq', 'ss'])
    comp_seq = str(seqres_data[base+"_"+key].seq)
    # Add this option for the cases where the IDX from the AA in PDB is not the
    # same of the original sequence
    seq_ss = ''.join(chain_df.loc[:,'seq'])
    new_seq_ss = get_alignment(comp_seq, seq_ss)
    non_gaps = np.array([x.start() for x in re.finditer('\w', new_seq_ss)])+1
    first_last_coords = str(non_gaps[0]-1)+"-"+str(chain_df.index[0])
    chain_df = chain_df.reset_index(drop=True).set_index(pd.Index(non_gaps))
    chain_df = chain_df.reindex(range(1, len(comp_seq)+1), fill_value="X")
    seq_mapped = ''.join(chain_df.loc[:,'seq'])
    ss_dssp = re.sub('-',' ',''.join(chain_df.loc[:,'ss']))
    new_lst.append({base+":"+key+":sequence:coords="+first_last_coords: seq_mapped})
    new_lst.append({base+":"+key+":secstr:coords="+first_last_coords: ss_dssp})
    return new_lst


def get_dssp_missing(comppath, base, seqres_data, un_pdb_ids):
    ext = os.path.splitext(os.path.basename(comppath))[1].replace(".", "")
    ss_new_lst, error_lst = list(), list()
    try:
        if (ext=="pdb"):
            pdb_parser = PDBParser(QUIET=True)
        else:
            pdb_parser = MMCIFParser()
        structure = pdb_parser.get_structure(base, comppath)
        model = structure[0]
        dssp = DSSP(model, comppath)
        ss_vals = np.empty([len(dssp),4], dtype=object)
        i=0
        for a_key in dssp.keys():
            ss_vals[i,0] = a_key[1][1]
            ss_vals[i,1] = a_key[0]
            ss_vals[i,2] = dssp[a_key][1]
            ss_vals[i,3] = dssp[a_key][2]
            i+=1
        ss_vals = ss_vals[np.lexsort((ss_vals[:,0], ss_vals[:,1]))]
        _ , ss_chain_start = np.unique(ss_vals[:,1], return_index=True)
        ss_chain_end = np.append(ss_chain_start[1:]-1, len(ss_vals)-1)
        for j in range(len(ss_chain_end)):
            chain_name = base+"_"+ss_vals[ss_chain_start[j], 1]
            if (un_pdb_ids is None) or (chain_name in un_pdb_ids): 
                chain_part = ss_vals[ss_chain_start[j]:ss_chain_end[j]+1, :]
                ss_new_lst = ss_new_lst+get_ss_string(chain_part, base, seqres_data)
            else:
                error_lst.append([chain_name, 'not a target'])
    except Exception as e:
        error_lst.append([base, type(e).__name__+ " - " + e.args[0]])
    return ss_new_lst, error_lst


def annotate_missing(path_selected, seqres_data, un_pdb_ids):
    ''' Generate the DSSP mask for the selected sequences. '''
    # Getting the complete path from the target files
    sel_files = select_target(path_selected, un_pdb_ids)
        
    ss_lst, error_lst = list(), list()
    temp_path = path_selected+'temp/'
    try:
        if not os.path.exists(temp_path):
            os.makedirs(temp_path)
    except OSError:
        print ("Creation of the directory %s failed" % temp_path)
    i=0
    start_time = time.time()
    for name in sel_files:
        basename = os.path.splitext(os.path.basename(name))[0]
        base = os.path.splitext(basename)[0]
        comppath = temp_path+basename
        try:
            with gzip.open(name, 'r') as f_in, open(comppath, 'wb') as f_out:
                shutil.copyfileobj(f_in, f_out)
        except Exception as e:
            error_lst.append([base, type(e).__name__+ " - " + e.args[0]])
        ss_new, error_val = get_dssp_missing(comppath, base, seqres_data, un_pdb_ids)
        error_lst = error_lst + error_val
        ss_lst = ss_lst + ss_new
        i+=1
        if (i%100==0):
            print("{0}% processed {1} seconds".format(str(round((i/len(sel_files)*100), 3)), str(time.time() - start_time)))
    try:
        shutil.rmtree(temp_path)
    except OSError as e:
        print("Error: %s : %s" % (temp_path, e.strerror))
    print("{0} seconds".format(time.time() - start_time))
    return ss_lst, error_lst



def main(comp_path):
    seqres_data, seqres_keys, seqres_desc, seqres_fnb = extract_seqres(path_seqres)
    missing_names, missing_ids, obsolete_names = get_missing_names(seqres_keys, path_masked)
    file_name = os.path.join(path_pdb_files, "missing_pdbs.txt")
    missing = manage_missing(path_pdb_files, missing_ids, file_name)
    if len(missing)>0:
        file_name = os.path.join(comp_path, datetime.today().strftime('%Y%m%d')+"_missing_pdbs_log.txt")
        resources.save_sep_file(missing, file_name)
        print("ATENTION: There are still some PDB/CIF files we were not able to download.\nCheck the file missing_pdbs_log and try to download them manually.\nWe will move on now considering the files available on disk!")
    ss2append, errors_dssp = annotate_missing(path_pdb_files, seqres_data, missing_ids)
