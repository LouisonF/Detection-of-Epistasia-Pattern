#! /usr/bin/python3
#-*-coding: utf-8-*-
# In order to launch the script in a conda environment, use this call: " Python ./file_to_launch"
# Louison Fresnais M2BB
# François Courtin M2BB

import os, os.path
import sys
import re

caus_snp_array = sys.argv[1].split(',')

TP = 0
FP = 0
FN = 0

list = os.listdir('../Genetic_results')
number_files = len(list)

for i in range(1, number_files+1):
    file  = open('../Genetic_results/' + str(i)+".txt", "r")
    #print(file.read())

    for line in file:
        inc_TP = False
        inc_FP = False
        inc_FN = False
        try:
            pattern = re.search('{(.+?)}', line).group(1)
        except AttributeError:
            pattern = ''
        if pattern == '':
            inc_FN = True
        else:
            snp_array = pattern.split(",")
            if (snp_array == caus_snp_array) or (snp_array[::-1] == caus_snp_array):
                inc_TP = True
            elif len(snp_array) < len(caus_snp_array):
                inc_FP = True
            else:
                for i in range(len(snp_array)):
                    for j in range(len(caus_snp_array)):
                        if snp_array[i] == caus_snp_array[j]:
                            inc_TP = True
            if not(inc_TP) and not(inc_FN):
                inc_FP = True


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

f_measure = 2*((TP/(FP+TP))*(TP/(TP+FN)))/((TP/(FP+TP))+(TP/(TP+FN)))

print("F-measure : ",f_measure)
