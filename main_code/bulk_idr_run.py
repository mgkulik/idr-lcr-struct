#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun  2 09:17:41 2023

@author: magoncal
"""

import os, time, gc
import idr, poly, idrAlphafold as idraf, alphafoldDssp as afdssp
import resources


comp_path = '/home/magoncal/Documents/data/projects/poly_cook/'
folders = sorted([f.name for f in os.scandir(comp_path) if os.path.isdir(f) and f.name.startswith('UP')])
#folders.remove('UP000005640')

start_time_all = time.time()

for un_prot in folders:
    
    comp_path_un = os.path.join(comp_path, un_prot)
    # MobiDB parameters
    path1 = os.path.join(comp_path_un, un_prot+'.mobi')
    key = "prediction-disorder-mobidb_lite"
    name = "mobidb"
    group="lite"
    # Fasta file
    path_fasta = sorted([f.name for f in os.scandir(comp_path_un) if f.name.endswith('.fasta')])
    #if len(path_fasta)==1:
    path_fasta = os.path.join(comp_path_un, path_fasta[0])
    
    # IDR output
    # start_time = time.time()
    tab_idr = os.path.join(comp_path_un, un_prot+"_mobidb_"+group+".tab")
    # idrs_path, idr_fasta_path = idr.run_all(path_fasta, tab_idr, path1, key, name, group)
    # print("Filtered fasta ({0}) and IDR details ({1}) were saved to disk.".format(os.path.basename(idr_fasta_path), os.path.basename(idrs_path)))
    # end_time = time.time()
    # time_formated = resources.transform_time(start_time, end_time)
    # print("{0} IDR FINISHED. Time: {1}".format(un_prot, time_formated))
    
    # # IDR properties output
    # start_time = time.time()
    # _ = idr.get_cider_props(idrs_path, "idr")
    # end_time = time.time()
    # time_formated = resources.transform_time(start_time, end_time)
    # print("{0} PROPS FINISHED. Time: {1}".format(un_prot, time_formated))
    
    # Poly output
    # start_time = time.time()
    
    print('\nStarting Poly steps:\n')
    # cutoff=.6
    # min_size=4
    # idrs_path = resources.gen_filename(tab_idr, "mobidb", group+"_details", "csv")
    
    # n_aa=1
    # source="polyx"
    # tab_poly = os.path.join(comp_path_un, un_prot+"_polyx.tab")
    # _ = poly.main_poly(path_fasta, tab_poly, idrs_path, source, int(n_aa), float(cutoff), int(min_size))
    
    # n_aa=2
    # source="polyxy"
    # tab_poly = os.path.join(comp_path_un, un_prot+"_polyxy.tab")
    # _ = poly.main_poly(path_fasta, tab_poly, idrs_path, source, int(n_aa), float(cutoff), int(min_size))
    
    # end_time = time.time()
    # time_formated = resources.transform_time(start_time, end_time)
    # print("{0} POLY FINISHED. Time: {1}".format(un_prot, time_formated))
    
    # AlphaFold step
    start_time = time.time()
    
    print('\nStarting AF steps:\n')
    afdssp.run_noselection(comp_path_un, un_prot, path_fasta)
    idraf.main_af_all(comp_path_un, un_prot)
    
    #else:
    #    print("\n{0} already generated".format(un_prot))
    
    
    end_time = time.time()
    time_formated = resources.transform_time(start_time, end_time)
    print("\n{0} AF FINISHED. Time: {1}".format(un_prot, time_formated))

end_time_all = time.time()
time_formated = resources.transform_time(start_time_all, end_time_all)
print("\nALL COMPLETE. Time: {0}".format(time_formated))