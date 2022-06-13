#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun  9 17:43:05 2022

This main file will receive different parameters according to what data
extraction should be performed.
TODO: Transforme this set of scripts on an object oriented package for 
future submission to PyPi.

@author: magoncal
"""

import idr
import os

seq = [1,5]
locals_var =  locals()


print("Please define what step of the process you wish to execute:")
print("1: Extract IDRs from MOBIDB;\n2: Generate IDR additional data;\n3: Generate PolyX/Ys additional data;\n4: Cross IDRs with PDB sequences;\n5: Cross Polys with PDB sequences.")
num_sel = input("Select a number according with description above ({0}-{1}): ".format(str(seq[0]), str(seq[1])))
assert(num_sel.isnumeric()), "Enter a number between {0} and {1}!".format(str(seq[0]), str(seq[1]))
assert(int(num_sel)>=seq[0] and int(num_sel)<=2), "Enter a number between {0} and {1}!".format(str(seq[0]), str(seq[1]))

if (int(num_sel)==1):
    """ Here your output will be a tab separated file with IDR regions and the 
    MobiDB-lite consensus score for the complete sequence. """
    
    print("Provide the complete path for the MobiDB json data.\nYou can download the organism complete proteome with: \nhttps://mobidb.bio.unipd.it/api/download?proteome=UP*****")
    path1 =  input("MobiDB json data: ")
    tab_file = idr.extract_json(path1)
    print("\nTab file {0} with MobiDB predictions was saved to disk.".format(os.path.basename(tab_file)))

elif (int(num_sel)==2):
    """ Here your output will be a csv file with several extractions and 
    calculations of IDR properties and a reduced fasta containing only the
    sequences with predicted IDRs. """
    
    path2 = input("Provide the complete organism fasta path: ")
    if not "tab_file" in locals_var:
        tab_file = input("Provide the IDR tab file path: ")
    
    idr_details, idr_fasta = idr.run_all(path2, tab_file)
    print("\nFiltered fasta ({0}) and IDR details ({1}) were saved to disk.".format(os.path.basename(idr_fasta), os.path.basename(idr_details)))
    
else:
    pass