#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun  9 17:43:05 2022

@author: magoncal
"""

"""
    This main file will receive different parameters according to what data
    extraction should be performed.
    TODO: Transforme this set of scripts on an object oriented package for 
    future submission to PyPi.
"""

seq = [1,4]


print("Please define what step of the process you wish to execute:")
print("1: Extract IDRs from MOBYDB;\n2: Annotate Polys;\n3: Cross IDRs with PDB sequences;\n4: Cross Polys with PDB sequences.")
num_sel = input("Select a number according with description above ({0}-{1}): ".format(str(seq[0]), str(seq[1])))
assert(num_sel.isnumeric()), "Enter a number between {0} and {1}!".format(str(seq[0]), str(seq[1]))
assert(int(num_sel)>=seq[0] and int(num_sel)<=2), "Enter a number between {0} and {1}!".format(str(seq[0]), str(seq[1]))

