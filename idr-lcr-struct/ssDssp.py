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

import resources

path_seqres = "/home/magoncal/Documents/data/projects/poly_cook/2022_07_11_pdb_seqres.txt"
path_pdb_files = "/home/magoncal/Documents/data/projects/idr_cook/pdb_files_missing_ss"
path_masked = "/home/magoncal/Documents/data/projects/poly_cook/2021_12_06_pdb_seqres_masked.txt"


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
        missing_ids = np.unique([n.split(sep)[0] for n in missing_names])
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


def manage_missing(path_pdb_files, missing_ids, file_name):
    ''' Deletes old PDB/CIF files that showed errors previously so the script
    can try to run DSSP again and download the .pdb file using the bash file 
    provided by PDB and when not available the .cif file. Every 10 executions 
    it checks if it should abort in case some files remain missing.
    '''
    
    # Removing files already existing in the folder (There may be errors from 
    # last execution).
    sel_files = select_target(path_pdb_files, list(missing_ids))
    if len(sel_files)>0 :
        try:
            for file in sel_files:
                os.remove(file)
        except Exception as e:
            print([base, e])
    
    # Now download all missing .pdb.gz files from PDB
    source_sh = os.path.join(path_pdb_files, "batch_download.sh")
    # Save missing files to disk, so the .sh provided by PDB can run over all 
    # downloads at once.
    file_name = os.path.join(path_pdb_files, "missing_pdbs.txt")
    missing = list(missing_ids)
    resources.save_sep_file(missing, file_name)
    prev_miss = len(missing)
    file_type = "-p"
    i = 0
    while len(missing) > 0:
        i += 1
        subprocess.call(shlex.split(f"{source_sh} -f {file_name} -o {path_pdb_files} {file_type}"))
        missing = list(select_target(path_pdb_files, missing, False))
        resources.save_sep_file(missing, file_name)
        if len(missing)!=prev_miss:
            prev_miss = len(missing)
        else:
            file_type = "-c"
        if (i%10==0 and file_type=="-c" and len(missing) > 0:
            print()
            print("{0} Attemps were already made and there are still missing PDB/CIF files.\nPlease download them manually.".format(str(i)))
            break
    os.remove(file_name)
    return missing
        
    
    


def main(comp_path):
    seqres_data, seqres_keys, seqres_desc, seqres_fnb = extract_seqres(path_seqres)
    missing_names, missing_ids, obsolete_names = get_missing_names(seqres_keys, path_masked)
    
