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

import os, subprocess, shlex, requests
from datetime import datetime
from gemmi import cif

import resources

#pdb_seq_url = "https://ftp.wwpdb.org/pub/pdb/derived_data/pdb_seqres.txt.gz"
#path_pdb_files = "/home/magoncal/Documents/data/projects/idr_cook/pdb_files_missing_ss"
#path_masked = "/home/magoncal/Documents/data/projects/poly_cook/pdb/20211206_pdb_seqres_masked.txt"
#path_ss_file = "/home/magoncal/Documents/data/projects/poly_cook/pdb/20211206_pdb_ss_final.txt"


def download_seqs_file(pdb_seq_url, path_seqres):
    url = pdb_seq_url
    response = requests.get(url)
    open(path_seqres+".gz", "wb").write(response.content)
    extract_from_gz(path_seqres+".gz", path_seqres)


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
    ''' Extract the list of files not available on the last extraction based on
    the masked fasta file.
    If no path_masked is provided we assume this is the first execution, so all
    protein PDB files will be downloaded.
    '''
    missing_names, missing_ids, obsolete_names = [], [], []
    
    if (path_masked!=""):
        ss_dict = resources.load_seqs(path_masked, "|", 'pdb')
        old_ss_keys = np.array(list(ss_dict.keys()))
        missing_names = np.setdiff1d(seqres_keys, old_ss_keys).tolist()
        missing_ids = np.unique([n.split(sep)[0] for n in missing_names]).tolist()
        obsolete_names = np.setdiff1d(old_ss_keys, seqres_keys).tolist()
    return missing_names, missing_ids, obsolete_names


def select_target(path_selected, un_pdb_ids, intersect=True, sep="_"):
    ''' Gets all files available on disk and extract just the complete path of 
    the target files based on an external list of IDs. 
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


def manage_missing(path_pdb_files, missing_ids, file_name, del_prev=True):
    ''' Deletes old PDB/CIF files that showed errors previously so the script
    can try to run DSSP again and download the .pdb file using the bash file 
    provided by PDB and when not available the .cif file. Every 10 executions 
    it checks if it should abort in case some files remain missing to avoid 
    infinite loops. It will usually finish in the third execution when the 
    internet connection is OK.
    '''
    
    # Removing files already existing in the folder (There may be errors from 
    # last execution). It will make the execution slower, but I found one case
    # where the PDB chain simply changed, but the ID was the same (6tif_AAA and 6tif_BBB)
    # became (6tif_A and 6tif_B), so I could have problems annotating again without a
    # new download.
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
    if len(missing) > 0:
        resources.save_sep_file(missing, file_name, "\n")
    prev_miss = len(missing)
    file_type = "-p"
    i = 0
    start_time = time.time()
    while len(missing) > 0:
        i += 1
        subprocess.call(shlex.split(f"{source_sh} -f {file_name} -o {path_pdb_files} {file_type}"))
        missing = list(select_target(path_pdb_files, missing, False))
        resources.save_sep_file(missing, file_name, "\n")
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
    if (i>0):
        os.remove(file_name)
    missing = list(select_target(path_pdb_files, missing_ids, False))
    return missing


def extract_from_gz(name, comppath):
    error = ""
    try:
        with gzip.open(name, 'r') as f_in, open(comppath, 'wb') as f_out:
            shutil.copyfileobj(f_in, f_out)
    except Exception as e:
        error = [base, type(e).__name__+ " - " + e.args[0]]
    return error
    

def get_one_pdb(base, ext, comppath):
    ''' Download again the PDB file. Needed when the download failed to save
    the correct pdb file to disk. '''
    # Deleting old files
    try:
        os.remove(comppath)
    except Exception:
        pass # If it can't delete it can move on.
    simpl_path = os.path.dirname(comppath).replace("/temp", "")
    simpl_path_file = os.path.join(simpl_path, base)+"."+ext+".gz"
    os.remove(simpl_path_file)
    
    # Downloading new one
    url = "https://files.rcsb.org/download/"+base+"."+ext+".gz"
    response = requests.get(url)
    open(simpl_path_file, "wb").write(response.content)
    
    # Extract from gz
    extract_from_gz(simpl_path_file, comppath)


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


def get_dssp_missing(comppath, base, seqres_data, un_pdb_names):
    ext = os.path.splitext(os.path.basename(comppath))[1].replace(".", "")
    ss_new_lst, error_lst = list(), list()
    eof=False
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
            if (len(un_pdb_names)==0 or (chain_name in un_pdb_names)): 
                chain_part = ss_vals[ss_chain_start[j]:ss_chain_end[j]+1, :]
                ss_new_lst = ss_new_lst+get_ss_string(chain_part, base, seqres_data)
            else:
                error_lst.append([chain_name, 'not a target'])
    except EOFError:
        eof=True
    except Exception as e:
        error_lst.append([base+"."+ext, type(e).__name__+ " - " + e.args[0]])
    return ss_new_lst, error_lst, eof


def annotate_missing(path_selected, seqres_data, un_pdb_ids, un_pdb_names=[]):
    ''' Generate the DSSP mask for the selected sequences. '''
    # Getting the complete path from the target files
    sel_files = select_target(path_selected, un_pdb_ids)
        
    ss_lst, error_lst = list(), list()
    temp_path = os.path.join(path_selected, 'temp')
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
        comppath = os.path.join(temp_path, basename)
        error_val = extract_from_gz(name, comppath)
        if len(error_val):
            error_lst.append(error_val)
        ss_new, error_val, eof = get_dssp_missing(comppath, base, seqres_data, un_pdb_names)
        if eof: # Deals with broken downloads
            get_one_pdb(base, ext, comppath)
            ss_new, error_val, eof = get_dssp_missing(comppath, base, seqres_data, un_pdb_names)
        error_lst = error_lst + error_val
        ss_lst = ss_lst + ss_new
        i+=1
        if (i%100==0):
            print()
            print("{0}% processed {1} seconds".format(str(round((i/len(sel_files)*100), 3)), str(time.time() - start_time)))
    try:
        shutil.rmtree(temp_path)
    except OSError as e:
        print("Error: %s : %s" % (temp_path, e.strerror))
    print("{0} seconds".format(time.time() - start_time))
    return ss_lst, error_lst


def read_fasta_ss(path_ss):
    temp = []
    with open (path_ss, "r") as file_handler:
        for line in file_handler:
            if line[:1] == ">":
                if temp:
                    comp_seq = re.sub('\\n', '', ''.join(temp))
                    temp = []
                    yield comp_seq, header
                header = line[1:]
            elif line[:1] != "\n":
                temp.append(line)
        if temp:
            comp_seq = re.sub('\\n', '', ''.join(temp))
            temp = []
            yield comp_seq, header


def extract_ss(path_ss, sep="_", s_type="pdb"):
    """Read the fasta file with the two dinstinct types of data (seq and SS)
       **PS: Decided not to use the SeqIO function because it deletes the
       white spaces from the dssp structures."""
    ss_data = dict()
    ss_seq = dict()
    ss_keys = list()
    ss_coords = dict()
    for v, k in read_fasta_ss(path_ss):
        ss_comp = k.split(':')
        ss_type = ss_comp[2]
        ss_name = ss_comp[0]
        if (s_type=="pdb"):
            ss_name+=sep+ss_comp[1]
        if (ss_type == 'secstr\n') or (ss_type == 'secstr'):
            ss_data[ss_name] = v
            ss_keys.append(ss_name)
            if len(ss_comp)>3:
                ss_coords[ss_name] = re.sub("coords=", "", ss_comp[3]).rstrip()
        else:
            ss_seq[ss_name] = v
    ss_keys = np.array(ss_keys)
    return ss_data, ss_keys, ss_seq, ss_coords
    
    
def remove_obsoletes(ss_obsolete, missing_ids, path_ss_file, path_new_ss, path_pdb_files):
    ''' Create a new multifasta with seqs and ss that should be kept,
    dealing with the obsolete entries from PDB. '''
    # Get data from old file to append to a new one
    ss_data, ss_keys, ss_seq, _ = extract_ss(path_ss_file)

    # select the ones to keep
    ss_keep = np.setdiff1d(ss_keys, ss_obsolete)
    
    # Add to the new file
    with open(path_new_ss, 'w') as newfile:
        for n in ss_keep:
            seq_name = re.sub('_', ':', n)
            newfile.write('>'+seq_name+':sequence\n')
            newfile.write(wrap_line(ss_seq[n],75))
            newfile.write('>'+seq_name+':secstr\n')
            newfile.write(wrap_line(ss_data[n],75))
    
    # Removing obsolete PDB files from disk
    pdb_obsolete = np.unique([ob.split("_")[0] for ob in ss_obsolete]).tolist()
    pdb_obsolete = np.setdiff1d(pdb_obsolete, missing_ids)
    for name in pdb_obsolete:
        path_file = os.path.join(path_pdb_files, name)
        if os.path.exists(path_file+".pdb.gz"):
            os.remove(path_file+".pdb.gz")
        if os.path.exists(path_file+".cif.gz"):
            os.remove(path_file+".cif.gz")


def wrap_line(seq, sz_wrap=75):
    ''' As I am saving the file directly without using biopyhton I need to 
    break the lines to generate the fasta.
    '''
    sz = len(seq)
    new_seq = '\n'.join([seq[i:i+sz_wrap] for i in range(0, sz, sz_wrap)])+'\n'
    return new_seq


def append_ss_tofasta(ss_fasta, path_new_ss):
    ''' Append the new sequences with secondary structures to the fasta file.
    IMPORTANT: THE RESULTING FILE IS UNORDERED.'''
    with open(path_new_ss, 'a+') as newfile:
        for l in ss_fasta:
            k = list(l.keys())[0]
            newfile.write('>'+ k +'\n')
            newfile.write(wrap_line(l[k],75))


def generate_masked_fasta(path_ss_final, path_masked, seqres_desc):
    ''' Extract the original PDB sequence with complete description to use as
    source of the blast database. '''
    struct2D = dict() 
    _, _, ss_seq, ss_coords = extract_ss(path_ss_final)
    with open(path_masked, 'w') as newfile:
        for k, seq in ss_seq.items():
            seq_name = re.sub('_', '|', k)
            newfile.write('>pdb|'+seq_name+' '+seqres_desc[k]+'\n')
            newfile.write(wrap_line(seq, 75))
            #coords=""
            #if k in ss_coords:
            #    coords=ss_coords[k]
            #struct2D[k] = [new_2Dstr, start_end, coords]
    #return struct2D
    


### Additional info from PDB files ###

def get_dbref_synthetic(comppath, base):
    ''' Direct reading of the PDB file structure to extract just the information
    of if the structure is syntetic and its global resolution. '''
    with open(comppath, 'r') as handle:
        rowsplit = ""
        base_synt = ""
        val_res = []
        for line in enumerate(handle):
            rowsplit = line[1].rstrip("\n").split(" ")
            rowsplit = [l for l in rowsplit if l!=""]
            # Find the synthetic proteins
            if rowsplit[0]=="SOURCE":
                content = ""
                idx_syntetic, idx_yes, idx_org_synt = -1, -1, -1
                content = line[1].rstrip("\n")
                idx_syntetic = content.find("SYNTHETIC")
                idx_yes = content.find("YES")
                idx_org_synt = content.find("SYNTHETIC CONSTRUCT")
                if idx_org_synt==-1:
                    idx_org_synt = content.find("ARTIFICIAL GENE")
                if idx_org_synt==-1:
                    idx_org_synt = content.find("32630")
                if (idx_syntetic!=-1 and idx_yes!=-1) or idx_org_synt!=-1:
                    base_synt = base
            if (rowsplit[0]=="REMARK"):
                if len(rowsplit)>2:
                    if (rowsplit[2]=="RESOLUTION."):
                        val = rowsplit[3]
                        if (val!="NOT" and val!="NULL"):
                            val_res.append([base, float(val)])
                if (rowsplit[1]=="3"):
                    break
    return base_synt, val_res


def get_dbref_synthetic_cif(comppath, base):
    ''' Extracts the information of if the structure is syntetic and its resolution
    from the cif file using the gemni package. '''
    base_synt = ""
    val = ""
    val_res = []
    with open(comppath, 'r') as handle:
        content = cif.read_string(handle.read())
    category_dict = content[0].get_mmcif_category('_pdbx_entity_src_syn')
    if (len(category_dict)!=0):
        key = "ncbi_taxonomy_id"
        for v in category_dict[key]:
            if (v is not None):
                ids = v.split(",")
                for i in ids:
                    if int(i)==32630:
                        base_synt = base
    if (base_synt==""):
        key = "pdbx_gene_src_ncbi_taxonomy_id"
        category_dict = content[0].get_mmcif_category('_entity_src_gen')
        if (len(category_dict)!=0):
            for v in category_dict[key]:
                if (v is not None):
                    ids = v.split(",")
                    for i in ids:
                        if int(i)==32630:
                            base_synt = base
    # Extracting the resolution from each aligned PDB
    category_dict = content[0].get_mmcif_category('_reflns')
    if (len(category_dict)!=0):
        val = category_dict["d_resolution_high"][0]
    if (isinstance(val, bool) or (val is None) or (val=="")):
        category_dict = content[0].get_mmcif_category('_refine')
        if (len(category_dict)!=0):
            val = category_dict["ls_d_res_high"][0]
    if (isinstance(val, bool) or (val is None) or (val=="")):
        category_dict = content[0].get_mmcif_category('_em_3d_reconstruction')
        if (len(category_dict)!=0):
            val = category_dict["resolution"][0]
    if ((val!="") or (val is not None)):
        val_res.append([base, float(val)])
    return base_synt, val_res


def extract_from_pdb_file_synt(path_pdb_files, sel_ids, resolution_path):
    ''' Function will extract all data from the PDB file in a first interaction
    then not used again for the future executions. Process added after the DSSP
    annotation. It would take to long to re-process all DSSP just to generate
    the first version of this file. '''
    start_time = time.time()
    ss_synt, ll_res =  [], []
    sel_files = np.array(sorted([f.path for f in os.scandir(path_pdb_files) if (re.search(r'\W*.gz', f.path))]))
    base_files = np.array([os.path.splitext(os.path.splitext(os.path.basename(f))[0])[0] for f in sel_files])
    bool_idx = np.isin(base_files, sel_ids)
    sel_files = sel_files[bool_idx].tolist()
    temp_path = path_pdb_files+'temp/'
    try:
        if not os.path.exists(temp_path):
            os.makedirs(temp_path)
    except OSError:
        print ("Creation of the directory %s failed" % temp_path)
    i=1
    for name in sel_files:
        out1, out2 = "", ""
        lst_jrnl, lst_res = [], []
        basename = os.path.splitext(os.path.basename(name))[0]
        base = os.path.splitext(basename)[0]
        ext = os.path.splitext(basename)[1]
        comppath = temp_path+basename
        with gzip.open(name, 'r') as f_in, open(comppath, 'wb') as f_out:
            shutil.copyfileobj(f_in, f_out)
        if ext==".pdb":
            out1, lst_res = get_dbref_synthetic(comppath, base)
        else:
            out1, lst_res = get_dbref_synthetic_cif(comppath, base)
        if (out1!=""):
            ss_synt.append(out1)
        if len(lst_res)>0:
            ll_res = ll_res + lst_res
        if (i%10**4==0):
            print("{0} of {1} completed for the synthetics search.".format(str(i), str(len(sel_files))))
        i+=1
    print("{0} completed for the synthetics search.".format(str(len(sel_files))))
    # Getting the resolution info
    df_res = pd.DataFrame(ll_res, columns=["pdb_id", "pdb_resolution", "syntetic"])
    tot_in_sec = time.time() - start_time
    print("--- %s seconds ---" % (tot_in_sec))
    try:
        shutil.rmtree(temp_path)
    except OSError as e:
        print("Error: %s : %s" % (temp_path, e.strerror))
    return ss_synt, df_res
    
    
def extract_from_pdb_file(path_pdb_files, path_pdb, date_start):
    ''' Searches for synthetics and experiment info for PDB IDs aligned to sequences . '''
    resolution_path = os.path.join(path_pdb, date_start+"_pdb_details.csv")
    if os.path.exists(resolution_path):
        df_resol = pd.read_csv(resolution_path)
    else:
        df_pdb = extract_from_pdb_file_synt(path_pdb_files, sel_ids, resolution_path)
        # This 4 PDB sequences have no organism and no indication of being synthetic.
        # They are old and engeneered. Checked even the .cif file, where they present an "?"
        lst_synt = lst_synt + ["1MJ0", "1N0Q", "1N0R", "3F0H"]
        save_list(ANALYSIS_PATH+synth_name, lst_synt)
        df_resol.to_csv(ANALYSIS_PATH+res_name, index=False)
    return lst_synt, df_resol


def main(path_pdb):
    
    date_start = datetime.today().strftime('%Y%m%d')
    path_seqres = os.path.join(path_pdb, date_start+"_pdb_seqres.txt")
    download_seqs_file(pdb_seq_url, path_seqres)
    seqres_data, seqres_keys, seqres_desc, seqres_fnb = extract_seqres(path_seqres)
    missing_names, missing_ids, obsolete_names = get_missing_names(seqres_keys, path_masked)
    
    file_name = os.path.join(path_pdb, date_start+"_missing_pdbs_list.txt")
    resources.save_sep_file(missing_ids, file_name)
    file_name = os.path.join(path_pdb_files, "missing_pdbs.txt")
    missing = manage_missing(path_pdb_files, missing_ids, file_name)
    if len(missing)>0:
        file_name = os.path.join(path_pdb, date_start+"_sitll_missing_pdbs_list.txt")
        resources.save_sep_file(missing, file_name)
        print("ATENTION: There are still some PDB/CIF files we were not able to download.\nCheck the file missing_pdbs_log and try to download them manually.\nWe will move on now considering the files available on disk!")
    
    # Dealing with the new file now. We remove the obsolete first and in the 
    # same function create the new file.
    path_ss_final = os.path.join(path_pdb, date_start+"_pdb_ss_final.txt")
    remove_obsoletes(obsolete_names, missing_ids, path_ss_file, path_ss_final, path_pdb_files)
    
    if len(missing_ids)>0:
        ss2append, errors_dssp = annotate_missing(path_pdb_files, seqres_data, missing_ids)
        # Save the transformation just to make sure, as it can take hours
        file_name = os.path.join(path_pdb, date_start+"_pdb_ss.pickle")
        resources.save_pickle(ss2append, file_name, 'wb')
        append_ss_tofasta(ss2append, path_ss_final)
        file_error = os.path.join(path_pdb, date_start+"_dssp_error.txt")
        resources.save_sep_llists(errors_dssp, file_error, "\t")
        path_masked_final = os.path.join(path_pdb, date_start+"_pdb_masked.txt")
        generate_masked_fasta(path_ss_final, path_masked_final, seqres_desc)
        #path_masked_pkl = os.path.join(path_pdb, date_start+"_pdb_masked.pickle")
        #resources.save_pickle(struct2D, path_masked_pkl, 'wb')
    
    # Extracting essential info from the PDB files (type, resolution, interactions?)
    # _, ss_keys, _, _ = extract_ss(path_ss_final)
    # sel_ids = np.unique([ids.split("_")[0] for ids in ss_keys]).tolist()  # Run this if you want to take all ids from 2D structure file
    df_pdb_details = extract_from_pdb_file(path_pdb_files, path_pdb, date_start)
   
    
    
