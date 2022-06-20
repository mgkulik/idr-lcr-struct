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
import poly

seq = [0,5]
locals_var =  locals()


print("Please define what step of the process you wish to execute:")
print("0: Execute all steps;\n1: Extract IDRs from MOBIDB;\n2: Generate IDR additional data;\n3: Generate PolyX/Ys additional data;\n4: Cross IDRs with PDB sequences;\n5: Cross Polys with PDB sequences.")
num_sel = input("Select a number according with description above ({0}-{1}): ".format(str(seq[0]), str(seq[1])))
assert(num_sel.isnumeric()), "Enter a number between {0} and {21}!".format(str(seq[0]), str(seq[1]))
assert(int(num_sel)>=seq[0] and int(num_sel)<=seq[1]), "Enter a number between {0} and {1}!".format(str(seq[0]), str(seq[1]))

if (int(num_sel)==1):
    """ Here your output will be a tab separated file with IDR regions and the 
    MobiDB-lite consensus score for the complete sequence. """
    
    print("All file names must start with the 11 characters uniprot proteome ID: e.g. UP000005640")
    print("Provide the complete path for the MobiDB json data.\nYou can download the organism complete proteome with: \nhttps://mobidb.bio.unipd.it/api/download?proteome=UP*****")
    path1 =  input("MobiDB json data: ")
    tab_idr = idr.extract_json(path1)
    print("\nTab file {0} with MobiDB predictions was saved to disk.".format(os.path.basename(tab_idr)))

elif (int(num_sel)==2):
    """ Here your output will be a csv file with several extractions and 
    calculations of IDR properties and a reduced fasta containing only the
    sequences with predicted IDRs. """
    
    path_fasta = input("Provide the path for the proteome fasta: ")
    if not "tab_idr" in locals_var:
        tab_idr = input("Provide the path for the IDR tab file: ")
    
    idr_details, idr_fasta = idr.run_all(path_fasta, tab_idr)
    print("\nFiltered fasta ({0}) and IDR details ({1}) were saved to disk.".format(os.path.basename(idr_fasta), os.path.basename(idr_details)))

elif (int(num_sel)==3):
    """ Here your output will be a csv file with several extractions and 
    calculations of PolyX or PolyXY properties. """
    
    cutoff=.6
    min_size=4
    
    tab_poly = input("Provide the path for the Poly tab file: ")   
    if not "path_fasta" in locals_var:
        path_fasta = input("Provide the path for the proteome fasta: ")
        
    if not "idr_details" in locals_var:
        idr_details = input("Provide the path for the IDR basic data (generated in step 2): ")
    
    n_aa = input("Number of different residues per repeat (e.g. 1=homorepeat, 2=direpeat ...): ")
    assert(n_aa.isnumeric()), "Value must be bigger than 0!"
    assert(int(n_aa)>0), "Value must be bigger than 0!"
    
    
    change_cut = input("The cut off values for poly inside IDR are: 0.60 or 4 residues.\nDo you want to change it (0:No, 1:Yes)? ")
    assert(change_cut.isnumeric()), "Value must be 0 or 1!"
    assert(int(change_cut)==0 or int(change_cut)==1), "Value must be 0 or 1!"
    
    if (int(change_cut)==1):
        cutoff = input("Define new accepted fraction (e.g. 0.60): ")
        assert(float(cutoff)>=0.5 and float(cutoff)<=1.0), "Value must be between 0 and 1!"
        min_size = input("Define new accepted minimum of residues (e.g. 4): ")
        assert(min_size.isnumeric()), "Value must be bigger than 4!"
        assert(int(min_size)>4), "Value must be bigger than 4!"
    
    poly_details = poly.run_poly(path_fasta, tab_poly, idr_details, int(n_aa), float(cutoff), int(min_size))
    print("IDR details ({0}) were saved to disk.".format(os.path.basename(poly_details)))

    
elif (int(num_sel)==4):
    
    if not "idr_details" in locals_var:
        idr_details = input("Provide the path for the IDR basic data (generated in step 2): ")
    
    
else:
    pass