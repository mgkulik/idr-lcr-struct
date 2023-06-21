#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 16 10:13:41 2022

@author: magoncal
"""

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.PDB import PDBParser, MMCIFParser
from Bio.PDB.DSSP import DSSP

from gemmi import cif

import pandas as pd
import numpy as np
import math
import re
import gzip, shutil, os, csv, pickle, shutil, time, tarfile

from datetime import datetime

import resources

#/home/magoncal/Documents/data/projects/idr_cook/
#basis_path = input("BE CAREFULL!! Set the complete path to the folder of the BASE files: ")
comp_path = '/home/magoncal/Documents/data/projects/poly_cook/'
un_prot = 'UP000005640'
comp_path_un = os.path.join(comp_path, un_prot)
path_fasta = os.path.join(comp_path_un, un_prot+'_9606.fasta')
#path_ids = basis_path+'uniprots_final_dataset'
#path_selected = '/home/magoncal/Documents/data/projects/poly_cook/UP000005640/'


def get_ss_string(chain_part, base, seqres_data, fold_id):
    ''' We are just assuming the sequence provided by alphafold initially because
    of differences caused by uniprot sequences changes. As there are no differences
    between their sequence and the final sequence, we will merge the long sequences
    later. '''
    new_lst = list()
    key = chain_part[0,1]
    chain_df = pd.DataFrame(chain_part[:,2:], index=chain_part[:,0], columns=['seq', 'ss'])
    # Add this option for the cases where the IDX from the AA in PDB is not the
    # same of the original sequence
    seq_mapped = ''.join(chain_df.loc[:,'seq'])
    ss_dssp = re.sub('-',' ',''.join(chain_df.loc[:,'ss']))
    new_lst.append({base+":"+fold_id+":"+"sequence": seq_mapped})
    new_lst.append({base+":"+fold_id+":"+"secstr": ss_dssp})
    return new_lst


def get_pLDDT_cif(comppath):
    with open(comppath, 'r') as handle:
        content = cif.read_string(handle.read())
    category_dict = content[0].get_mmcif_category('_atom_site')
    if (len(category_dict)!=0):
        key = "auth_seq_id"
        id_pos = np.array(category_dict[key])
        key = "B_iso_or_equiv"
        val_pos = np.array(category_dict[key])
        all_vals = np.vstack((id_pos, val_pos))
        _, idxs = np.unique(all_vals, axis=1, return_index=True)
        all_vals = all_vals[:, np.sort(idxs)]
    return (all_vals[1,:])


def get_dssp_missing(comppath, base, seqres_data, file_name, fold_id):
    ext = os.path.splitext(os.path.basename(comppath))[1].replace(".", "")
    ss_new_lst, error_lst = list(), list()
    try:
        if (ext=="pdb"):
            pdb_parser = PDBParser(QUIET=True)
        else:
            pdb_parser = MMCIFParser()
        structure = pdb_parser.get_structure(file_name, comppath)
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
            chain_part = ss_vals[ss_chain_start[j]:ss_chain_end[j]+1, :]
            ss_new_lst = ss_new_lst+get_ss_string(chain_part, base, seqres_data, fold_id)
    except Exception as e:
        error_lst.append([base, type(e).__name__+ " - " + str(e.args[0])])
    return ss_new_lst, error_lst


def sort_files(sel_files):
    # Extract the model numbers
    sel_files = np.sort(sel_files)
    sel_models = np.array([os.path.splitext(os.path.splitext(os.path.basename(f))[0])[0].split("-")[2].replace("F", "") for f in sel_files], dtype="int64")
    # Generate a numeric representation of the uniprot ID to use as group in the sorting
    base_files = np.array([os.path.splitext(os.path.splitext(os.path.basename(f))[0])[0].split("-")[1] for f in sel_files])
    _, base_idx = np.unique(base_files, return_inverse=True)
    base_rep = np.arange(0, len(base_idx))[base_idx]
    base_rep = np.column_stack((base_rep, sel_models))
    # Extract the ordered index and apply in the original array
    idx_sort = np.lexsort((base_rep[:, 1], base_rep[:, 0]))
    sel_files = sel_files[idx_sort].tolist()
    #sel_names = np.array([os.path.splitext(os.path.splitext(os.path.basename(f))[0])[0].split("-") for f in sel_files])
    return sel_files


def select_target(path_selected, un_pdb_ids, ptype=""):
    files = sorted([f.path for f in os.scandir(path_selected) if (re.search(r'\W*'+ptype+'.gz', f.path))])
    base_files = [os.path.splitext(os.path.splitext(os.path.basename(f))[0])[0].split("-")[1] for f in files]
    un_pdb = np.unique([n.split("_")[0] for n in un_pdb_ids])
    sel_files = []
    for x, y in zip(base_files, files):
        if (x in un_pdb_ids):
            sel_files.append(y)
    sel_files = sort_files(np.array(sel_files))
    return(sel_files)


def gen_str_row(base, plddt):
    ''' Merge the ID and values in a tab file. '''
    return(base + '\t' + ', '.join(list(plddt)))


def extract_selected(comp_path_un, fasta_selected, ptype="cif"):
    error_lst = []
    tar_path = [(f.path, f.name) for f in os.scandir(comp_path_un) if f.name.endswith('.tar')][0]
    # It always create a folder inside another. Managing this
    cif_path = tar_path[0][:-4]
    cif_subpath = os.path.join(cif_path, tar_path[1][:-4])
    
    files_tar = tarfile.open(tar_path[0], 'r')
    for prot_name in fasta_selected.keys():
        miss = True
        for member in files_tar.getmembers():
            if prot_name in member.name:
                files_tar.extract(member, cif_path)
                miss = False
                break
        if miss:
            error_lst.append([prot_name, 'Not available in the AlphaFold tar file.'])
    return cif_subpath, cif_path, os.listdir(cif_subpath), error_lst


def extract_from_tar(comp_path_un, ptype="cif"):
    tar_path = [(f.path, f.name) for f in os.scandir(comp_path_un) if f.name.endswith('.tar')][0]
    # It always create a folder inside another. Managing this
    cif_path = tar_path[0][:-4]
    cif_subpath = os.path.join(cif_path, tar_path[1][:-4])
    files_tar = tarfile.open(tar_path[0], 'r')
    for member in files_tar.getmembers():
        files_tar.extract(member, cif_path)
    return cif_subpath, cif_path, np.array(os.listdir(cif_subpath))


def annotate_ss(comp_path_un, fasta_selected, ptype="cif", dssp=True):
    ss_lst, plddt_lst, error_lst = list(), list(), list()
    # Dropped the idea of selecting the files because the search is too much time consuming
    #cif_subpath, cif_path, sel_files, error_lst = extract_selected(comp_path_un, fasta_selected, ptype)
    # Decided to filter just the ones present in the fastar later
    cif_subpath, cif_path, sel_files = extract_from_tar(comp_path_un, ptype)
    base_files = [n.split("-")[1] for n in sel_files]
    bool_idx = np.isin(base_files, list(fasta_selected.keys()))
    sel_files = sel_files[bool_idx]
    sel_files = [os.path.join(cif_subpath, f) for f in sel_files]
    # Sorting list based on the number of the model (now it is considering numbers as char)
    sel_files = sort_files(sel_files)

    i=0
    start_time = time.time()
    for name in sel_files:
        basename = os.path.splitext(os.path.basename(name))[0]
        file_name = os.path.splitext(basename)[0]
        base = file_name.split("-")[1]
        fold_id = file_name.split("-")[2].replace("F", "")
        comppath = os.path.join(cif_subpath, basename)
        try:
            with gzip.open(name, 'r') as f_in, open(comppath, 'wb') as f_out:
                shutil.copyfileobj(f_in, f_out)
        except Exception as e:
            error_lst.append([base, type(e).__name__+ " - " + str(e.args[0])])
        if dssp:
            ss_new, error_val = get_dssp_missing(comppath, base, fasta_selected, file_name, fold_id)
            error_lst = error_lst + error_val
            ss_lst = ss_lst + ss_new
        plddt = get_pLDDT_cif(comppath)
        plddt_lst.append(gen_str_row(base+"_"+fold_id, plddt))
        os.remove(comppath)
        i+=1
        if (i%1000==0):
            print("{0}% processed {1} seconds".format(str(round((i/len(sel_files)*100), 3)), str(time.time() - start_time)))
    try:
        pass
        shutil.rmtree(cif_path)
    except OSError as e:
        print("Error: %s : %s" % (temp_path, e.strerror))
    print("{0} seconds".format(time.time() - start_time))
    return ss_lst, plddt_lst, error_lst


def wrap_line(seq, sz_wrap=75):
    sz = len(seq)
    new_seq = '\n'.join([seq[i:i+sz_wrap] for i in range(0, sz, sz_wrap)])+'\n'
    return new_seq


def append_ss_tofasta(ss_fasta, path_ss, name):
    with open(path_ss, 'a+') as newfile:
        for l in ss_fasta:
            k = list(l.keys())[0]
            newfile.write('>'+ k +'\n')
            newfile.write(wrap_line(l[k],75))

            
def select_fasta(path_fasta):
    ''' Reads the fasta file and store the sequences.
    Returns a SeqIO object to use later. '''
    
    fasta_lst = []
    fasta_data = dict()
    for seq_record in SeqIO.parse(path_fasta, 'fasta'):
        uniprot_name = seq_record.id.split('|')[1]
        fasta_lst.append(seq_record)
        fasta_data[uniprot_name] = seq_record
    
    return fasta_data


def save_list(path, lst_synt):
    textfile = open(path, "w")
    for element in lst_synt:
        textfile.write(element + "\n")
    textfile.close()


def run_noselection():
    start_time = time.time()
    # I decided to annotate all the proteome and filter the sequences with IDRs later
    fasta_selected = select_fasta(path_fasta)
    #fasta_selected_ori = fasta_selected.copy(); fasta_selected = {k: fasta_selected[k] for k in list(fasta_selected)[:10]}
    print("Starting DSSP step ...\n")
    ss2append, plddts, errors_dssp = annotate_ss(comp_path_un, fasta_selected, "cif") #Almost 5hs
    if len(errors_dssp)>0:
        file_error = os.path.join(comp_path_un, un_prot+"_alphafold_error.txt")
        resources.save_sep_llists(errors_dssp, file_error, "\t")
    path_ss = os.path.join(comp_path_un, un_prot+'_ss_alphafold.fasta')
    append_ss_tofasta(ss2append, path_ss, "new")
    save_list(os.path.join(comp_path_un, un_prot+"_plddts.tab"), plddts)
    with open(os.path.join(comp_path_un, un_prot+'_ss_alphafold.pickle'), 'wb') as handle:
        pickle.dump(ss2append, handle, protocol=pickle.HIGHEST_PROTOCOL)
    end_time = time.time()
    time_formated = resources.transform_time(start_time, end_time)
    print("\nPART 1 FINISHED. Time: {0}".format(time_formated))