#! /usr/bin/python3
#-*-coding: utf-8-*-
# In order to launch the script in a conda environment, use this call: " Python ./file_to_launch"
# Louison Fresnais M2BB
# François Courtin M2BB

#run : ./set_gamet.py input_directory_path output_directory_path

import os, os.path
import sys

input = sys.argv[1]
output = sys.argv[2]


for gamete in os.listdir(input):
    pheno = []
    geno = ""
    file = open(input+"/"+gamete, "r")
    geno_file = open(output+"/geno_"+os.path.basename(gamete), "a")

    for line in file:
        ar_line  = line.split('\t')
        pheno.append(ar_line[len(ar_line)-1])
        ar_line.pop(len(ar_line)-1)
        geno_file.write(','.join(ar_line)+'\n')

    pheno_file = open(output+"/pheno_"+os.path.basename(gamete), "a")
    pheno_file.write(''.join(pheno))
