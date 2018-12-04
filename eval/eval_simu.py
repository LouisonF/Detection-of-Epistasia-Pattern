#! /usr/bin/python3
#-*-coding: utf-8-*-
# In order to launch the script in a conda environment, use this call: " Python ./file_to_launch"
# Louison Fresnais M2BB
# François Courtin M2BB

#run : ./f_measure input_directory_path output_directory_path n_runs

from utils import *
import os, os.path
import sys
import re
import glob






"""
    if (inc_TP):
        TP += 1

    elif(inc_FN):
        FN += 1
    elif(inc_FP):
        FP += 1

print("TP : ",TP)
print("FP : ",FP)
print("FN : ",FN)

if (FP == 0 and FN != 0):
    FP = 1
if (FN == 0 and FP != 0):
    FN = 1

#recall = TP/(TP+FN)
#precision = TP/(TP+FP)
#f_measure = 2/((1/recall)+(1/precision))

#print("F-measure : ",f_measure)
"""
