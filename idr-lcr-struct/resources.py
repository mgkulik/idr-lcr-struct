#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 13 16:20:12 2022

Collection of auxiliary functions used by several source codes.

@author: magoncal
"""
import os
import pickle
from Bio import SeqIO
import numpy as np

def save_file(lst, path):
    ''' Save list to txt file '''
    with open(path, 'a') as file:
        for data in lst:
          file.write(data)
          
def gen_filename(filename, source, proc, ext):
    ''' Extracts the info from the source file name and generates a new file name.'''
    dir_name = os.path.dirname(filename)
    file_name = os.path.basename(filename).split('.')[0].split('_')[0]
    new_name = dir_name+"/"+file_name+"_"+source+"_"+proc+"."+ext
    return new_name

def save_pickle(var, var_name, BASIS_PATH):
    with open(BASIS_PATH+var_name, 'ab') as handle:
        pickle.dump(var, handle, protocol=pickle.HIGHEST_PROTOCOL)

def open_pickle(var_name, BASIS_PATH):
    objs = []
    with open(BASIS_PATH+var_name, 'rb') as handle:
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


def read_fasta(fastaname, seq_ints=False):
    ''' Reads the fasta file and store the sequences and their int equivalents.
    Counts AA by groups and total of residues of all proteome. '''
    sz = get_filesize(fastaname, ">")
    seqs_by_group = np.zeros((sz,4))
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


def save_fastas(seq_lst_comp, fastaout):
    ''' Saving filtered fasta files to use in other tasks. '''    
    with open(fastaout, 'w') as handle:
      SeqIO.write(seq_lst_comp, handle, "fasta")