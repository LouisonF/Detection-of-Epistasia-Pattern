#! /usr/bin/python3
#-*-coding: utf-8-*-
# In order to launch the script in a conda environment, use this call: " Python ./file_to_launch"
# Louison Fresnais M2BB
# François Courtin M2BB

#run : ./f_measure input_directory_path output_directory_path n_runs causal_snp

from utils import *
import os, os.path
import sys
import re

input = sys.argv[1]
output = sys.argv[2]
n_runs = sys.argv[3]
caus_snp = sys.argv[4]


#List of directories of conditionnal dataset
for condition in os.listdir(input):

    os.mkdir(output+"/"+condition)
    print(os.path.basename(condition))
    #List of directories of simulated conditionnal dataset
    for simu_file in os.listdir(input+"/"+condition):

        #List of n_runs result files for one simulated conditionnal dataset
        for res_file in os.listdir(input+"/"+condition+"/"+simu_file):
            #Get TP or FP or FN for each result file
            result = set_res(input+"/"+condition+"/"+simu_file+"/"+res_file, caus_snp)

            #Write the result in the f_measure.txt file
            output_file = write_res(output+"/"+condition+"/"+simu_file, result, "result")

        print(os.path.basename(output_file))
        f_measure = run_f_measure(output_file)
        write_res(output+"/"+condition, str(f_measure), "f_measure")

        power = run_power(output_file, int(n_runs))
        write_res(output+"/"+condition, str(power), "power")
