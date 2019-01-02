#! /usr/bin/python3
#-*-coding: utf-8-*-
# In order to launch the script in a conda environment, use this call: " Python ./file_to_launch"
# Louison Fresnais M2BB
# François Courtin M2BB

#run : ./eval_simu.py input_directory_path output_directory_path n_runs nb_snp len_pattern

from utils import *
import os, os.path
import sys
import re

input = sys.argv[1]
output = sys.argv[2]
n_runs = sys.argv[3]
nb_snp = sys.argv[4] #number of snp in the matrix
len_pattern = sys.argv[5]
caus_snp = ""

for i in range(int(len_pattern)):
    caus_snp += "N"+str(int(nb_snp) - int(i))+","

caus_snp = caus_snp[0:(len(caus_snp)-1)] #pop the "," at the end of the string

print("causal SNPs" + caus_snp)
#List of directories of conditionnal dataset

os.mkdir(output+"/"+os.path.basename(input))
condition = os.path.basename(output+"/"+os.path.basename(input))
print(os.path.basename(input))
#List of directories of simulated conditionnal dataset
for simu_file in os.listdir(input):

    #List of n_runs result files for one simulated conditionnal dataset
    for res_file in os.listdir(input+"/"+simu_file):
        #Get TP or FP or FN for each result file
        result = set_res(input+"/"+simu_file+"/"+res_file, caus_snp)

        #Write the result in the f_measure.txt file
        output_file = write_res(output+"/"+condition+"/"+simu_file, result, "result")

    print(os.path.basename(output_file))
    f_measure = run_f_measure(output_file)
    write_res(output+"/"+condition, str(f_measure), "f_measure")

    power = run_power(output_file, int(n_runs))
    write_res(output+"/"+condition, str(power), "power")
