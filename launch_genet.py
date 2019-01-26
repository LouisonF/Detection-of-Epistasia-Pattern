#! /usr/bin/python3
#-*-coding: utf-8-*-
# In order to launch the script in a conda environment, use this call: " Python ./file_to_launch"
# Louison Fresnais M2BB
# François Courtin M2BB

import os, os.path
import sys

geno_path = sys.argv[2]
pheno_path = sys.argv[3]

for geno in os.listdir(geno_path):
    for i in range(1,21):
        pheno = "pheno" + geno[4:]
        print(geno)
        print(pheno)
        os.system("./Genetic/Debug/Genetic res_"+os.path.basename(geno)+"_"+str(i)+" "+geno_path+"/"+geno+" "+pheno_path+"/"+pheno)
