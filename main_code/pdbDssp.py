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
from gemmi import cif
import pandas as pd
import numpy as np
import math
import re


import gzip, shutil, pickle, time
import os, subprocess, shlex, requests
from gzip import BadGzipFile
from datetime import datetime


import resources

#pdb_seq_url = "https://ftp.wwpdb.org/pub/pdb/derived_data/pdb_seqres.txt.gz"
#path_pdb_files = "/home/magoncal/Documents/data/projects/idr_cook/pdb_files_missing_ss"
#path_masked = "/home/magoncal/Documents/data/projects/poly_cook/pdb/20211206_pdb_seqres_masked.txt"
#path_ss_file = "/home/magoncal/Documents/data/projects/poly_cook/pdb/20211206_pdb_ss_final.txt"


def download_seqs_file(pdb_seq_url, path_seqres):
    ''' Downloads the PDB active sequences file directly from PDB. '''
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


def select_target(path_selected, un_pdb_ids, file_type_name="", intersect=True, sep="_"):
    ''' Gets all files available on disk and extract just the complete path of 
    the target files based on an external list of IDs. 
            un_pdb_ids must have the PDB ID WHITOUT CHAIN.
            intersect: Searches for the intersection or the difference between sets?
    '''
    pattern = r'\W*'+file_type_name+'.gz'
    files = sorted([f.path for f in os.scandir(path_selected) if (re.search(pattern, f.path))])
    base_files = [os.path.splitext(os.path.splitext(os.path.basename(f))[0])[0] for f in files]
    # Remove duplicates when both cif and pdb are in the file, keeping the .cif
    _, check_dup = np.unique(base_files, return_index=True)
    files = [files[i] for i in check_dup]
    base_files = [base_files[i] for i in check_dup]
    
    sel_files = []
    if (intersect):
        for x, y in zip(base_files, files):
            if (x in un_pdb_ids):
                sel_files.append(y)
    else:
        sel_files = np.setdiff1d(un_pdb_ids, base_files)
    return sel_files


def manage_missing(path_pdb_files, missing_ids, step="-p", del_prev=True):
    ''' Deletes old PDB/CIF files that showed errors previously so the script
    can try to run DSSP again and download the .pdb file using the bash file 
    provided by PDB and when not available the .cif file. Every 10 executions 
    it checks if it should abort in case some files remain missing to avoid 
    infinite loops. It will usually finish in the third execution when the 
    internet connection is OK.
    
        step: Define what kind of files will be dowloaded as first option of
        pdb files (-p or -c);
        del_prev: Should previous files be deleted from the folder?
    '''
    file_type_name = ""
    if step=="-c":
        file_type_name = ".cif"
        
    file_name = os.path.join(path_pdb_files, "missing_pdbs.txt")
    
    # Removing files already existing in the folder (There may be errors from 
    # last execution). It will make the execution slower, but I found one case
    # where the PDB chain simply changed, but the ID was the same (6tif_AAA and 6tif_BBB)
    # became (6tif_A and 6tif_B), so I could have problems annotating again without a
    # new download.
    if del_prev:
        sel_files = select_target(path_pdb_files, missing_ids, file_type_name)
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
    missing = list(select_target(path_pdb_files, missing_ids, file_type_name, False))
    if len(missing) > 0:
        resources.save_sep_file(missing, file_name, "\n")
    prev_miss = len(missing)
    file_type = step
    i = 0
    start_time = time.time()
    while len(missing) > 0:
        i += 1
        subprocess.call(shlex.split(f"{source_sh} -f {file_name} -o {path_pdb_files} {file_type}"))
        missing = list(select_target(path_pdb_files, missing, file_type_name, False))
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
    missing = list(select_target(path_pdb_files, missing_ids, file_type_name, False))
    return missing


def extract_from_gz(name, comppath):
    ''' Extract txt file with gzip. '''
    error = ""
    try:
        with gzip.open(name, 'r') as f_in, open(comppath, 'wb') as f_out:
            shutil.copyfileobj(f_in, f_out)
    except Exception as e:
        error = [type(e).__name__+ " - " + e.args[0]]
    return error
    

def get_one_pdb(base, ext, comppath):
    ''' Download again the PDB file. Needed when the download failed to save
    the correct pdb file to disk. '''
    # Deleting old files
    try:
        os.remove(comppath)
    except Exception:
        pass # If it can't delete it can move on.
    simpl_path = os.path.dirname(comppath).replace("temp", "")
    simpl_path_file = os.path.join(simpl_path, base)+ext+".gz"
    os.remove(simpl_path_file)
    
    # Downloading new one
    url = "https://files.rcsb.org/download/"+base+ext
    response = requests.get(url)
    gzip.open(simpl_path_file, 'wb').write(response.content)
    
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
        ext = os.path.splitext(basename)[1]
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
        os.remove(comppath)
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
  
    
def remove_obsoletes(ss_obsolete, missing_ids, path_ss_file, path_new_ss, path_pdb_files):
    ''' Create a new multifasta with seqs and ss that should be kept,
    dealing with the obsolete entries from PDB. '''
    # Get data from old file to append to a new one
    ss_data, ss_keys, ss_seq, _ = resources.extract_ss(path_ss_file)

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
    ss_data, _, ss_seq, ss_coords = resources.extract_ss(path_ss_final)
    with open(path_masked, 'w') as newfile:
        for k, seq in ss_seq.items():
            seq_name = re.sub('_', '|', k)
            newfile.write('>pdb|'+seq_name+' '+seqres_desc[k]+'\n')
            newfile.write(wrap_line(seq, 75))
            coords=""
            if k in ss_coords:
                coords=ss_coords[k]
            struct2D[k] = [ss_data[k], coords]
    return struct2D
    


### Additional info from PDB files ###

def get_dbref_synthetic(comppath, base):
    ''' Direct reading of the PDB file structure to extract just the information
    of if the structure is syntetic and its global resolution. '''
    with open(comppath, 'r') as handle:
        rowsplit = ""
        base_synt = ""
        rel_date = []
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
            if (rowsplit[0]=="REVDAT"):
                if len(rowsplit)>=5:
                    if (rowsplit[4]=="0"):
                        rel_date.append([base, datetime.strptime(rowsplit[2],'%d-%b-%y').strftime('%Y-%m-%d')])
            if (rowsplit[0]=="REMARK"):
                if len(rowsplit)>2:
                    if (rowsplit[2]=="RESOLUTION."):
                        val = rowsplit[3]
                        if (val!="NOT" and val!="NULL"):
                            val_res.append([base, float(val)])
                if (rowsplit[1]=="3"):
                    break
    return base_synt, val_res, rel_date


def get_dbref_synthetic_cif(comppath, base):
    ''' Extracts the information of if the structure is syntetic and its resolution
    from the cif file using the gemni package. '''
    base_synt = ""
    val = " "  # Had to add a space here because empty can be extracted from the data
    rel_date = []
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
    # Release date
    category_dict = content[0].get_mmcif_category('_pdbx_audit_revision_history')
    if (len(category_dict)!=0):
        val = category_dict["data_content_type"]
        if (val[0] == "Structure model"):
            rel_date.append([base, category_dict["revision_date"][0]])
        val = " "
    # Extracting the resolution from each aligned PDB
    category_dict = content[0].get_mmcif_category('_reflns')
    if (len(category_dict)!=0):
        val = category_dict["d_resolution_high"][0]
    if (isinstance(val, bool) or (val is None) or (val==" ")):
        category_dict = content[0].get_mmcif_category('_refine')
        if (len(category_dict)!=0):
            val = category_dict["ls_d_res_high"][0]
    if (isinstance(val, bool) or (val is None) or (val==" ")):
        category_dict = content[0].get_mmcif_category('_em_3d_reconstruction')
        if (len(category_dict)!=0):
            val = category_dict["resolution"][0]
    if ((val!=" ") and (val is not None) and (val is not None)):
        #print(base)
        val_res.append([base, float(val)])
    return base_synt, val_res, rel_date


def append_synt_df(df_res, ss_synt):
    ''' Takes the syntetic IDs and add to the dataframe with the resolution data. '''
    # This step was required because there are cases of PDB/CIF files without
    # resolution, so I needed to make sure any synthetic is added independently
    # of having a resolution or not.
    df_synt = pd.DataFrame(list(zip(ss_synt,[True]*len(ss_synt))), columns=['pdb_id','synthetic'])
    df_res = pd.merge(df_res, df_synt, how="outer", on="pdb_id")
    df_res["synthetic"].fillna(False,inplace=True)
    return df_res


def filter_pdb_files(path_pdb_files, sel_ids):
    ''' As Biopython DSSP annotations still have some limitations for cif files, 
    to save disk space we prioritize PDB files downloads over cif, but after
    the target selection we decided to use CIF files to extract extra info, because
    they are more reliable. This caused some duplications of IDs on disk, so we 
    must select our targets prioritizing the PDB files for general executions. '''
    # Get PDB files
    sel_files_pdb = np.array(sorted([f.path for f in os.scandir(path_pdb_files) if (re.search(r'\W*pdb.gz', f.path))]))
    base_files_pdb = np.array([os.path.splitext(os.path.splitext(os.path.basename(f))[0])[0] for f in sel_files_pdb])
    bool_idx_pdb = np.isin(base_files_pdb, sel_ids)
    # Get CIF files and select just the difference
    sel_files_cif = np.array(sorted([f.path for f in os.scandir(path_pdb_files) if (re.search(r'\W*cif.gz', f.path))]))
    base_files_cif = np.array([os.path.splitext(os.path.splitext(os.path.basename(f))[0])[0] for f in sel_files_cif])
    diff_cif = np.setdiff1d(base_files_cif, base_files_pdb)
    # Select only the target IDs
    dif_idx = np.isin(base_files_cif, diff_cif)
    sel_files_cif = sel_files_cif[dif_idx]
    base_files_cif = base_files_cif[dif_idx]
    bool_idx_cif = np.isin(base_files_cif, sel_ids)
    # Merge PDB and cif files
    #sel_names = base_files_pdb[bool_idx_pdb].tolist()
    #sel_names = sel_names+base_files_cif[bool_idx_cif].tolist()
    #sel_names.sort()==sel_ids.sort()
    sel_files = sel_files_pdb[bool_idx_pdb].tolist()
    sel_files = sel_files+sel_files_cif[bool_idx_cif].tolist()
    return (sel_files)


def extract_from_pdb_file_synt(path_pdb_files, sel_ids, resolution_path):
    ''' Function will extract all data from the PDB file in a first interaction
    then not used again for the future executions. Process added after the DSSP
    annotation. It would take to long to re-process all DSSP just to generate
    the first version of this file. '''
    start_time = time.time()
    ss_synt, ll_res, ll_release =  [], [], []
    sel_files = filter_pdb_files(path_pdb_files, sel_ids)
    temp_path = os.path.join(path_pdb_files, 'temp')
    try:
        if not os.path.exists(temp_path):
            os.makedirs(temp_path)
    except OSError:
        print ("Creation of the directory %s failed" % temp_path)
    i=1
    for name in sel_files:
        out1, sel_files = "", ""
        lst_res = []
        basename = os.path.splitext(os.path.basename(name))[0]
        base = os.path.splitext(basename)[0]
        ext = os.path.splitext(basename)[1]
        comppath = os.path.join(temp_path, basename)
        with gzip.open(name, 'r') as f_in, open(comppath, 'wb') as f_out:
            shutil.copyfileobj(f_in, f_out)
        if ext==".pdb":
            out1, lst_res, rel_date = get_dbref_synthetic(comppath, base)
        else:
            out1, lst_res, rel_date = get_dbref_synthetic_cif(comppath, base)
        if (out1!=""):
            ss_synt.append(out1)
        if len(lst_res)>0:
            ll_res = ll_res + lst_res
        if len(rel_date)>0:
            ll_release = ll_release+rel_date
        if (i%10**4==0):
            print("{0} of {1} completed for the synthetics search.".format(str(i), str(len(sel_files))))
        os.remove(comppath)
        i+=1
    print("{0} completed for the synthetics search.".format(str(len(sel_files))))
    # This 4 PDB sequences have no organism and no indication of being synthetic.
    # They are old and engeneered. Checked even the .cif file, where they present an "?"
    ss_synt = ss_synt + ["1MJ0", "1N0Q", "1N0R", "3F0H"]
    # Getting the resolution info and the release date
    df_res = pd.DataFrame(ll_res, columns=["pdb_id", "pdb_resolution"])
    df_rel_date = pd.DataFrame(ll_release, columns=["pdb_id", "pdb_release"])
    df_res = pd.merge(df_res, df_rel_date, how="outer", on="pdb_id")
    # Appending synthetic IDs
    df_res = append_synt_df(df_res, ss_synt)
    tot_in_sec = time.time() - start_time
    print("--- %s seconds ---" % (tot_in_sec))
    try:
        shutil.rmtree(temp_path)
    except OSError as e:
        print("Error: %s : %s" % (temp_path, e.strerror))
    return df_res


def extract_from_pdb_file(path_pdb_files, path_pdb, date_start, sel_ids, det_old_path=""):
    ''' Load or annotate synthetics and other experiment info for PDB IDs 
    aligned to sequences . '''
    resolution_path = os.path.join(path_pdb, date_start+"_pdb_details.csv")
    df_pdb_new = extract_from_pdb_file_synt(path_pdb_files, sel_ids, resolution_path)
    if os.path.exists(det_old_path):
        df_pdb = pd.read_csv(det_old_path)
        df_pdb = pd.concat([df_pdb_new, df_pdb]).sort_values(by=["pdb_id"])
    df_pdb.to_csv(resolution_path, index=False)
    return df_pdb


def extract_dbref_cif(comppath, base):
    ''' Open each cif file and extract the targeted information:
    Source for the PDB files:
        
    Equivalences between PDB and CIF files:
    https://mmcif.wwpdb.org/docs/pdb_to_pdbx_correspondences.html
    '''
    lst_db_ref, lst_chain, lst_organism, lst_exp = [], [], [], []
    eof=False
    try:
        with open(comppath, 'r') as handle:
            content = cif.read_string(handle.read())
        # Equivalent to DBREF, DBREF1 and DBREF2 (PDB coords)
        category_dict = content[0].get_mmcif_category('_struct_ref_seq')
        if (len(category_dict)!=0):
            for i in range(0, len(category_dict["align_id"])):
                pdb_id, pdb_ref_db, pdb_ref = "","",""
                pdb_start, pdb_end, pdb_start_ref, pdb_end_ref = 0,0,0,0
                pdb_id = base
                pdb_name = pdb_id+'_'+category_dict["pdbx_strand_id"][i]
                pdb_ref = category_dict["pdbx_db_accession"][i]
                if ("pdbx_auth_seq_align_beg" in category_dict):
                    if (category_dict["pdbx_auth_seq_align_beg"][i] != None) and (category_dict["pdbx_auth_seq_align_end"][i] != None):
                        pdb_start_auth = category_dict["pdbx_auth_seq_align_beg"][i]
                        pdb_end_auth = category_dict["pdbx_auth_seq_align_end"][i]
                    else:
                        pdb_start_auth = np.nan
                        pdb_end_auth = np.nan
                if (category_dict["seq_align_beg"][i] != None) and (category_dict["seq_align_end"][i] != None):
                    pdb_start = category_dict["seq_align_beg"][i]
                    pdb_end = category_dict["seq_align_end"][i]
                pdb_ref_db = "UNP"
                if (len(pdb_ref)==4):
                    pdb_ref_db = "PDB"
                if (re.match('^[1-9]\d*$', pdb_start) is None) or (re.match('^[1-9]\d*$', pdb_end) is None):
                    special_coords = True
                else:
                    special_coords = False
                pdb_start_ref = int(category_dict["db_align_beg"][i])
                pdb_end_ref = int(category_dict["db_align_end"][i])
                lst_db_ref.append([pdb_name, pdb_id, pdb_start, pdb_end, pdb_start_auth, pdb_end_auth,
                                   pdb_ref_db, pdb_ref, pdb_start_ref, pdb_end_ref, special_coords])
        # Equivalent to COMPND, to find the strand of the molecule
        # Decided to extract directly from the source because the documentation
        # describes it as 100% informed, although not mandatory
        category_dict = content[0].get_mmcif_category('_entity_poly')
        if (len(category_dict)!=0):
            for i in range(0, len(category_dict["entity_id"])):
                molid = "MOL_ID: "+category_dict["entity_id"][i]
                chains = category_dict["pdbx_strand_id"][i].replace(" ", "").split(",")
                lst_chain = lst_chain + np.transpose(np.vstack([[base]*len(chains), [molid]*len(chains), chains])).tolist()
        # Equivalent to SOURCE, to extract the organism ID
        # I was using first _gen and if it didn't exist _nat, but sometimes they happend together
        category_dict = content[0].get_mmcif_category('_entity_src_gen')
        if (len(category_dict)!=0):
            keys = ["pdbx_gene_src_ncbi_taxonomy_id", "pdbx_gene_src_scientific_name", "pdbx_beg_seq_num", "pdbx_end_seq_num"]
            for i in range(0, len(category_dict["entity_id"])):
                molid = "MOL_ID: "+category_dict["entity_id"][i]
                organism_id = category_dict[keys[0]][i]
                organism_name = category_dict[keys[1]][i]
                start_coords = category_dict[keys[2]][i]
                end_coords = category_dict[keys[3]][i]
                if (organism_name is not None) and (organism_id is not None):
                    lst_organism.append([base, molid, organism_id, organism_name.upper(), start_coords, end_coords, "gen"])
        category_dict = content[0].get_mmcif_category('_entity_src_nat')
        if (len(category_dict)!=0):
            keys = ["pdbx_ncbi_taxonomy_id", "pdbx_organism_scientific", "pdbx_beg_seq_num", "pdbx_end_seq_num"]
            for i in range(0, len(category_dict["entity_id"])):
                molid = "MOL_ID: "+category_dict["entity_id"][i]
                organism_id = category_dict[keys[0]][i]
                organism_name = category_dict[keys[1]][i]
                start_coords = category_dict[keys[2]][i]
                end_coords = category_dict[keys[3]][i]
                if (organism_name is not None) and (organism_id is not None):
                    lst_organism.append([base, molid, organism_id, organism_name.upper(), start_coords, end_coords, "nat"])
        # Equivalent to EXPDTA
        category_dict = content[0].get_mmcif_category('_exptl')
        if (len(category_dict)!=0):
            lst_exp.append([base, category_dict["method"][0]])
    except EOFError:
        eof=True
    except ValueError:
        # Experiment 7KW7 don't have a cif available. Ignoring and just showing the error
        eof=True
    return lst_db_ref, lst_chain, lst_organism, lst_exp, eof


def join_organisms_coords(df_organism, df_pdb_coords):
    ''' Found out that the chimeric sequences can belong to several different
    organisms. There is no direct link between organisms and aligned id, so 
    we need to use the provided coordinates to map them together.
    '''
    return df_pdb_coords

def extract_from_cif_all(path_pdb_files, sel_ids):
    ''' Get all extra information from the CIF files, now targeting just the 
    selected structures. Decided to use only CIF files because the AUTH field
    is always there, while the PDB file may or may not have the data. '''
    start_time = time.time()
    print("\nStarting the extraction of extra information from the cif files...")
    sel_files, ss_lst, ss_chain, ss_organism, ss_exp = [], [], [], [], []
    files = sorted([f.path for f in os.scandir(path_pdb_files) if (re.search(r'\W*.cif.gz', f.path))])
    base_files = [os.path.splitext(os.path.splitext(os.path.basename(f))[0])[0] for f in files]
    
    for x, y in zip(base_files, files):
        if (x in sel_ids):
            sel_files.append(y)
    temp_path = os.path.join(path_pdb_files, 'temp')
    try:
        if not os.path.exists(temp_path):
            os.makedirs(temp_path)
    except OSError:
        print ("Creation of the directory %s failed" % temp_path)
    i=0
    for name in sel_files:
        i += 1
        basename = os.path.splitext(os.path.basename(name))[0]
        base = os.path.splitext(basename)[0]
        ext = os.path.splitext(basename)[1]
        comppath = os.path.join(temp_path, basename)
        try:
            with gzip.open(name, 'r') as f_in, open(comppath, 'wb') as f_out:
                shutil.copyfileobj(f_in, f_out)
        except Exception:
            get_one_pdb(base, ext, comppath)
        ss_coords, lst_chain, lst_organism, lst_exp, eof = extract_dbref_cif(comppath, base)
        if eof:
            get_one_pdb(base, ext, comppath)
            ss_coords, lst_chain, lst_organism, lst_exp, eof = extract_dbref_cif(comppath, base)
        ss_lst = ss_lst + ss_coords
        ss_chain = ss_chain + lst_chain
        ss_organism = ss_organism + lst_organism
        ss_exp = ss_exp+lst_exp
        os.remove(comppath)
        if eof:
            print('{0}.cif.gz had a download problem.\nPlease download it again and place in the right folder.'.format(base))
            break
        if (i%1000==0):
            print("Data extracted from {0} of {1} CIF files.".format(str(i), str(len(sel_files))))
    # The differences can be generated for several reasons. I decided to use
    # the data just when the mapping do not start in the position 1 in DBREF
    cols=['pdb_name', 'pdb_id', 'dbrefs_start', 'dbrefs_end', 
          'dbrefs_auth_start', 'dbrefs_auth_end', 'pdb_ref_db', 'align_id', 
          'align_start', 'align_end', 'special_coords']
    df_pdb_coords = pd.DataFrame(ss_lst, columns=cols)
    df_chain = pd.DataFrame(ss_chain, columns=["pdb_id", "mol_id", "chain"])
    # Experiment types
    df_exp = pd.DataFrame(ss_exp, columns=["pdb_id", "experiment_desc"])
    # All files should have coordinates and experiments. Merging and not discarding any side
    df_pdb_coords = pd.merge(df_pdb_coords, df_exp, how="outer", on="pdb_id")
    # aligned organism data
    cols_org = ["pdb_id", "mol_id", "aligned_organism_id", "aligned_organism_name",
                "dbrefs_start", "dbrefs_end", "src"]
    df_organism = pd.DataFrame(ss_organism, columns=cols_org)
    # Need to use the coordinates provided in the organism session to find the
    # equivalent aligned ID (start and final coords can be different)
    df_pdb_coords = join_organisms_coords(df_organism, df_pdb_coords)
    tot_in_sec = time.time() - start_time
    print("\n--- %s seconds ---" % (tot_in_sec))
    try:
        shutil.rmtree(temp_path)
    except OSError as e:
        print("Error: %s : %s" % (temp_path, e.strerror))
    return df_pdb_coords, df_chain


def main_files(path_pdb):
    
    date_start = datetime.today().strftime('%Y%m%d')
    path_seqres = os.path.join(path_pdb, date_start+"_pdb_seqres.txt")
    download_seqs_file(pdb_seq_url, path_seqres)
    seqres_data, seqres_keys, seqres_desc, seqres_fnb = extract_seqres(path_seqres)
    missing_names, missing_ids, obsolete_names = get_missing_names(seqres_keys, path_masked)
    
    file_name = os.path.join(path_pdb, date_start+"_missing_pdbs_list.txt")
    resources.save_sep_file(missing_ids, file_name)
    missing = manage_missing(path_pdb_files, missing_ids)
    if len(missing)>0:
        file_name = os.path.join(path_pdb, date_start+"_still_missing_pdbs_list.txt")
        resources.save_sep_file(missing, file_name)
        print("ATENTION: There are still some PDB/CIF files we were not able to download.\nCheck the file _sitll_missing_pdbs_list.txt and try to download them manually.\nWe will move on now considering the files available on disk!")


def main_dssp(path_pdb, path_ss_file, path_pdb_files, path_seqres):
    
    seqres_data, seqres_keys, seqres_desc, seqres_fnb = extract_seqres(path_seqres)
    missing_names, missing_ids, obsolete_names = get_missing_names(seqres_keys, path_masked)
    
    date_start = os.path.basename(path_seqres).split("_")[0]
    
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
        struct2D = generate_masked_fasta(path_ss_final, path_masked_final, seqres_desc)
        path_masked_pkl = os.path.join(path_pdb, date_start+"_ss_masked.pickle")
        resources.save_pickle(struct2D, path_masked_pkl, 'wb')
    
        # Extracting essential info from the PDB files (type, resolution, interactions?)
        # _, ss_keys, _, _ = resources.extract_ss(path_ss_final)
        # sel_ids = np.unique([ids.split("_")[0] for ids in ss_keys]).tolist()  # Run this if you want to take all ids from 2D structure file
        # sel_ids = missing_ids.copy()
        # det_old_path = ""
        #df_pdb_details = extract_from_pdb_file(path_pdb_files, path_pdb, date_start, sel_ids, det_old_path)
        
        
