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
import idrPdb
import resources

seq = [0,6]
file_names=['evalues_idrs', 'idr_sizes', 'unique_ids', 'gen_info', 'blast_over']
locals_var =  locals()

if not "un_prot" in locals_var:
    un_prot = input("Please provide the uniprot proteome ID starting in UP and followed by 9 numbers: ")
    assert(resources.check_uniprot_name(un_prot)), "Provide an uniprot proteome ID in the format: UN00XXXXXXXX (2 letters and 9 numbers)."

if not "comp_path" in locals_var:
    comp_path = input("Inform the folder where all your files will be saved: ")
    assert(os.path.exists(comp_path)), "Make sure the script can save files in: {0}.".format(comp_path)
    assert(os.path.isdir(comp_path)), "Provide a directory, not a file!"
    
    # Create the organism directory to store the output (It gets a mess when several executions are made)
    comp_path_un = os.path.join(comp_path, un_prot)
    if not os.path.exists(comp_path_un):
        os.mkdir(comp_path_un)
        print("Please insert all the required files on the organism folder: {0}".format(comp_path_un))

print("Please define what step of the process you wish to execute:")
print("0: Generate masked PDB sequnces and 2D; \
      \n1: Execute all steps; \
      \n2: Extract IDRs from MOBIDB; \
      \n3: Generate IDR additional data; \
      \n4: Generate PolyX/Ys additional data; \
      \n5: Cross IDRs with PDB sequences; \
      \n6: Cross Polys with PDB sequences.")
num_sel = input("Select a number according with description above ({0}-{1}): ".format(str(seq[0]), str(seq[1])))
assert(num_sel.isnumeric()), "Enter a number between {0} and {21}!".format(str(seq[0]), str(seq[1]))
assert(int(num_sel)>=seq[0] and int(num_sel)<=seq[1]), "Enter a number between {0} and {1}!".format(str(seq[0]), str(seq[1]))


if int(num_sel)==0:
    ''' Download the newest list of PDB sequences, PDB/CIF files, annotate
    2D structures with DSSP and generate the masked file of PDB sequences based 
    on the missing residues. '''
    
    if not pdb_seq_url in locals_var:
        pdb_seq_url = "https://ftp.wwpdb.org/pub/pdb/derived_data/pdb_seqres.txt.gz"
        check = input("The current URL to download the PDB sequences is: {0}. Do you want to change it (y/n)?".format(pdb_seq_url))
        assert(str(check).upper()!="Y" or str(check).upper()!="N"), "y for yes and n for no!"
        if check.upper()=="Y":
            pdb_seq_url = input("New URL: ")
    
    if not "path_pdb_files" in locals_var:
        path_pdb_files = input("Please inform the complete path where the PDB/CIF files will be stored.\nIMPORTANT: Active PDB structure files require at least 50GB of disk space: ")
        assert(os.path.exists(path_pdb_files)), "Please provide an existing directory with the proper permissions to save and delete files."
        if path_pdb_files[-1]!="/":
            path_pdb_files = path_pdb_files+"/"
        
    pdb_files = resources.get_pdb_directory(comp_path)
    
    if not "path_masked" in locals_var:
        message = "No previous masked file located in the pdb folder.\nPlease confirm it, otherwise all PDB files will be downloaded.\\(your process will run for several days)."
        path_masked = resources.get_pdbOut_names(files, '_pdb_masked.txt', message)
       
    if not "path_ss_file" in locals_var:
        message = "No previous 2D structures file located in the pdb folder.\nPlease confirm it, otherwise all PDB files will be downloaded.\\(your process will run for several days)."
        path_ss_file = resources.get_pdbOut_names(files, '_pdb_ss_final.txt', message)


if (int(num_sel)==1 or int(num_sel)==2):
    """ Here your output will be a tab separated file with IDR regions and the 
    MobiDB-lite consensus score for the complete sequence. """
    
    print("All file names must start with the 11 characters uniprot proteome ID: e.g. UP000005640")
    print("Provide the filename for MobiDB json data.\nYou can download the organism complete proteome with: \nhttps://mobidb.bio.unipd.it/api/download?proteome=UP*****")
    if not "path1" in locals_var:
        path1 = input("Filename for MobiDB json data.\nIt must be available in the directory provided before and start with the Uniprot proteome ID: ")
        path1 = resources.valid_file(path1, comp_path_un)
    tab_idr = idr.extract_json(path1)
    print("\nTab file {0} with MobiDB predictions was saved to disk.".format(os.path.basename(tab_idr)))

if (int(num_sel)==1 or int(num_sel)==3):
    """ Here your output will be a csv file with several extractions and 
    calculations of IDR properties and a reduced fasta containing only the
    sequences with predicted IDRs. """
    
    if not "path_fasta" in locals_var:
        path_fasta = input("Provide the path for the proteome fasta. \nIt must be available in the directory provided before and start with the Uniprot proteome ID: ")
        path_fasta = resources.valid_file(path_fasta, comp_path_un)
    if not "tab_idr" in locals_var:
        tab_idr = os.path.join(comp_path_un, un_prot+"_mobidb_idr.tab")
    
    idr_details, idr_fasta = idr.run_all(path_fasta, tab_idr)
    print("\nFiltered fasta ({0}) and IDR details ({1}) were saved to disk.".format(os.path.basename(idr_fasta), os.path.basename(idr_details)))


if (int(num_sel)==1 or int(num_sel)==4):
    """ Here your output will be a csv file with several extractions and 
    calculations of PolyX or PolyXY properties. """
    
    cutoff=.6
    min_size=4
    
    if not "tab_poly" in locals_var:
        tab_poly = input("Provide the path for the Poly tab file. \nIt must be available in the directory provided before and start with the Uniprot proteome ID:  ")   
        tab_poly = resources.valid_file(tab_poly, comp_path_un)
    if not "path_fasta" in locals_var:
        path_fasta = input("Provide the path for the proteome fasta. \nIt must be available in the directory provided before and start with the Uniprot proteome ID: ")
        path_fasta = resources.valid_file(path_fasta, comp_path_un)
        
    if not "idr_details" in locals_var:
        idr_details = os.path.join(comp_path_un, un_prot+"_mobidb_idr_details.csv")
    
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
        assert(min_size.isnumeric()), "Value must be higher than 4!"
        assert(int(min_size)>4), "Value must be higher than 4!"
    
    poly_details = poly.run_poly(path_fasta, tab_poly, idr_details, int(n_aa), float(cutoff), int(min_size))
    print("IDR details ({0}) were saved to disk.".format(os.path.basename(poly_details)))

    
if (int(num_sel)==1 or int(num_sel)==5):
    """ Several outputs that save the crossing data do disc to save memory (using pickle files) .
    
    Main output will be a big csv file with all the cases where and IDR
    in our set overlapped with a PDB sequence in a significant way. All IDRs not
    overlapped with PDB are also kept for future analysis.
    IMPORTANT: This file may have multiple entries for each IDR. """
    
    cutoff_idr=.5
    min_size_idr=10
    
    if not "idr_details" in locals_var:
        idr_details = os.path.join(comp_path_un, un_prot+"_mobidb_idr_details.csv")
     
    if not "blast_path" in locals_var:
        blast_path = input("Provide the path for the blast XML file. \nIt must be available in the directory provided before and start with the Uniprot proteome ID: ")   
        blast_path = resources.valid_file(blast_path, comp_path_un)
    
    change_cut = input("The cut off values for IDRs overlapping PDB sequences are: 0.50 or 10 residues.\nDo you want to change it (0:No, 1:Yes)? ")
    assert(change_cut.isnumeric()), "Value must be 0 or 1!"
    assert(int(change_cut)==0 or int(change_cut)==1), "Value must be 0 or 1!"
    
    if (int(change_cut)==1):
        cutoff_idr = input("Define new accepted fraction (e.g. 0.50): ")
        assert(float(cutoff_idr)>=0.5 and float(cutoff_idr)<=1.0), "Value must be between 0 and 1!"
        min_size_idr = input("Define new accepted minimum of residues (e.g. 10): ")
        assert(min_size_idr.isnumeric()), "Value must be higher than 10!"
        assert(int(min_size_idr)>10), "Value must be higher than 10!"
    
    idrs_path = idrPdb.merge_idr_pdb(idr_details, blast_path, file_names, (cutoff_idr, min_size_idr))
    
    pdb_files = resources.get_pdb_directory(comp_path)
    
    if not "path_masked" in locals_var:
        message = "No previous masked file located in the pdb folder.\nYou can't execute this step before execute step 0."
        path_masked = resources.get_pdbOut_names(pdb_files, '_pdb_masked.txt', message)
    
    if not "path_ss_file" in locals_var:
        message = "No previous 2D structures file located in the pdb folder.\nYou can't execute this step before execute step 0."
        path_ss_file = resources.get_pdbOut_names(pdb_files, '_pdb_ss_final.txt', message)
        
    if not "dssp_path" in locals_var:
        message = "No previous 2D pickle file located in the pdb folder.\nYou can't execute this step before execute step 0."
        dssp_path = resources.get_pdbOut_names(pdb_files, '_pdb_ss.pickle', message)
    
    idrPdb.run_ss_annotation(idrs_path, path_masked, path_ss_file, dssp_path, save_dup=True)