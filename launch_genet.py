#! /usr/bin/python3
#-*-coding: utf-8-*-
# In order to launch the script in a conda environment, use this call: " Python ./file_to_launch"
# Louison Fresnais M2BB
# François Courtin M2BB

import os, os.path
import sys

geno_path = sys.argv[1]
pheno_path = sys.argv[2]
dataset = sys.argv[3] #nom du jeu de données
os.mkdir("Genetic_results/"+dataset)

for geno in os.listdir(geno_path):
    geno_dir = "Genetic_results/"+ dataset + "/" + geno[:-13]
    os.mkdir(geno_dir)
    for i in range(1,21):
        pheno = geno[:-12] + "Phenotype.txt"
        print(geno)
        print(pheno)
        os.system("./Genetic/Debug/Genetic "+geno_dir+"/res_"+os.path.basename(geno)+"_"+str(i)+" "+geno_path+"/"+geno+" "+pheno_path+"/"+pheno)
