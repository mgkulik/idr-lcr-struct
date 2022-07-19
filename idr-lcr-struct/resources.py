#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 13 16:20:12 2022

Collection of auxiliary functions used by several source codes.

@author: magoncal
"""
import os
import pickle
import csv
from Bio import SeqIO
import numpy as np
import pandas as pd

AA_CODE_LIST = ['?','A','R','N','D','C','Q','E','G','H','I','L','K','M','F','P','S','T','W','Y','V','U','O','B','Z','X']
AA_CODE_DICT = {'A':1,'C':2,'D':3,'E':4,'F':5,'G':6,'H':7,'I':8,'K':9,'L':10,'M':11,'N':12,'P':13,'Q':14,'R':15,'S':16,'T':17,'V':18,'W':19,'Y': 20,'U':21,'O':22,'B':23,'Z':24,'X':25}
AA_CODE_ABRV = ['?', 'Ala', 'Arg', 'Asn', 'Asp', 'Cys', 'Gln', 'Glu', 'Gly', 'His', 'Ile', 'Leu', 'Lys', 'Met', 'Phe', 'Pro', 'Ser', 'Thr', 'Trp', 'Try', 'Val', 'Sec', 'Pyl', 'Asx', 'Glx', 'Xaa']
AA_GROUPS = [[19,8,16,17,5,3,6,22], [14,18,1,15,20,11,10,13,21], [12,2,9], [4,7]]
AA_GROUPS_NAMES = ['Polar Uncharged', 'Non-Polar', 'Polar Basic', 'Polar Acidic']

# I decided to create this simplified polarity groups to use for poly and some other analysis
AA_GROUPS2 = [[19,8,16,17,5,3,6,22,12,2,9,4,7], [14,18,1,15,20,11,10,13,21]]
AA_GROUPS2_NAMES = ['Polar', 'Non-Polar']

def check_uniprot_name(name):
    ''' Simple name check to make sure a proteome name was provided and that
    all files start with the proteome name. '''
    name = name.split(".")[0].split("_")[0]
    return name[0:2] and name[2:].isnumeric() and len(name)==11


def valid_file(filename, path):
    comp_path = os.path.join(path, filename)
    assert(check_uniprot_name(filename)), "All file names must start with the uniprot proteome ID: UN00XXXXXXXX (2 letters and 9 numbers)."
    assert(os.path.isfile(comp_path)), "Provide a filename of an existent file in {0}!".format(path)
    return comp_path


def save_file(lst, path):
    ''' Save list to txt file with no separator. '''
    with open(path, 'a') as file:
        for data in lst:
          file.write(data)

def save_sep_file(lst, path, nrow="", sep=","):
    ''' Save a file with delimiter using csv structure from a simple list. 
    PS: I needed to add a row delimiter rdel option to the end of the file because
    PDB batch_download.sh doesn't with the last ID. '''
    with open(path, 'w', newline="") as myfile:
        wr = csv.writer(myfile, delimiter=sep)
        wr.writerow(lst+[nrow])


def save_sep_llists(lst, path, sep=","):
    ''' Save a file with delimiter using csv structure from a list of lists. '''
    with open(path, 'w', newline="") as myfile:
        wr = csv.writer(myfile, delimiter="\t")
        wr.writerows(lst)


def get_dir(filename):
    return os.path.dirname(filename)


def get_filename(filename):
    return os.path.basename(filename).split('.')[0].split('_')[0]


def gen_filename(filename, source, proc, ext):
    ''' Extracts the info from the source file name and generates a new file name.'''
    dir_name = get_dir(filename)
    file_name = get_filename(filename)
    sep = "_"
    if proc=="":
        sep = proc
    new_name = dir_name+"/"+file_name+"_"+source+sep+proc+"."+ext
    return new_name


def save_pickle(var, var_name, mode):
    with open(var_name, mode) as handle:
        pickle.dump(var, handle, protocol=pickle.HIGHEST_PROTOCOL)

def open_pickle(var_name, BASIS_PATH):
    objs = []
    with open(BASIS_PATH+var_name+".pickle", 'rb') as handle:
        while 1:
            try:
                objs.append(pickle.load(handle))
            except EOFError:
                break
    var = np.concatenate(objs, axis=0)
    return var

def get_filesize(file, item):
    ''' Counts the number of entries of a file to use to create numpy objects. '''
    return len([1 for line in open(file) if line.startswith(item)])


def get_seq_ints(seq):
    ''' Gets sequence AAs and generate their equivalent numeric sequence'''
    seq_int = np.array([AA_CODE_LIST.index(aa) for aa in seq])
    seq_by_aa = np.bincount(seq_int)
    seq_by_aa = np.pad(seq_by_aa, (0, len(AA_CODE_LIST) - len(seq_by_aa)), 'constant')
    seq_by_group = [sum(seq_by_aa[groups]) for groups in AA_GROUPS]
    return seq_by_group


def read_fasta(fastaname, seq_ints=False):
    ''' Reads the fasta file and store the sequences and their int equivalents.
    Counts AA by groups and total of residues of all proteome. 
    
    Returns a SeqIO object to use later. '''
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


def load_seqs(path, sep="|", s_type='seq'):
    '''Loads the sequence from the fasta file and returns a dictionary. '''
    seq_dict = dict()
    for seq_record in SeqIO.parse(path, 'fasta'):
        name_seq = seq_record.id.split(sep)
        if (s_type=='seq'):
            name_seq = name_seq[1]
        elif (s_type=='pdb'):
            #name_seq = name_seq[1].upper()+'_'+name_seq[2].upper()
            name_seq = name_seq[1]+'_'+name_seq[2]
        seq_dict[name_seq] = [str(seq_record.seq), seq_record.description]
    return seq_dict


def save_fastas(seq_lst_comp, fastaout):
    ''' Saving filtered fasta files to use in other tasks. '''    
    with open(fastaout, 'w') as handle:
      SeqIO.write(seq_lst_comp, handle, "fasta")


def extract_proteome(filename):
    ''' Extract code of the proteome from file name. '''
    return int(os.path.basename(filename).split('.')[0].split('_')[0][2:])


def get_aa_type(aa):
    ''' Get a simplified informaiton of the AA groups. '''
    int_aa = AA_CODE_LIST.index(aa)
    idx_type = [ind for ind, x in enumerate(AA_GROUPS) if int_aa in x][0]
    return AA_GROUPS_NAMES[idx_type]


def get_AAcomposition(seq_aa, tp=1, gp=AA_GROUPS2):
    ''' Counts the composition of AAs by group and returns lists or series
    according to how it is used (1:pandas or 2:list) and group  '''
    idr_int = [AA_CODE_LIST.index(aa) for aa in seq_aa]
    idr_group_tots = [sum([el in groups for el in idr_int]) for groups in gp]
    idr_group_props = [i/len(seq_aa) for i in idr_group_tots]
    if (tp==1):
        return pd.Series((v for v in idr_group_tots))
    else:
        return idr_group_tots, idr_group_props
    
                
def get_idr_AAcomposition(idr_aa):
    ''' Calculate the Physical-Chemical proportions of the IDR and extract first 
    and second most common group. '''
    idr_group_tots, idr_group_props = get_AAcomposition(idr_aa, 2, AA_GROUPS)
    idr_tops_diff = sorted(idr_group_tots, reverse=True)
    idr_tops_diff = idr_tops_diff[0]-idr_tops_diff[1]
    idr_tops_diff_prop = sorted(idr_group_props, reverse=True)
    idr_tops_diff_prop = idr_tops_diff_prop[0]-idr_tops_diff_prop[1]
    # Gets the variance between the groups
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


def get_pairs(col1, col2):
    ''' Reorder the poly pair identification to group the polyXYs composed of
    the same amino acids, e.g. "EK" == "KE"
    '''
    if AA_CODE_DICT[col1]<AA_CODE_DICT[col2]:
        val = col1+'|'+col2
    else:
        val = col2+'|'+col1
    return val


def extract_ss(path_ss, sep="_", s_type="pdb"):
    """Read the fasta file special for the secondary structure cases.
    Here we have 1 entry with the sequence (sequence) and 1 with the secondary
    structure labels (secstr).
       **PS: Decided not to use the SeqIO function because it deletes the
       white spaces from the dssp structures. """
    ss_data = dict()
    ss_seq = dict()
    ss_keys = list()
    ss_coords = dict()
    for v, k in read_fasta(path_ss):
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