#! /usr/bin/python3
#-*-coding: utf-8-*-
# launch example: python launch_smmb.py ../Simu_naive/Simu_naive_2snp_0.5_geno/ ../Simu_naive/Simu_naive_2snp_0.5_pheno/ Simu_naive_2snp_0.5
# Louison Fresnais M2BB
# François Courtin M2BB

import os, os.path
import sys

geno_path = sys.argv[1] #Give the genotypes directory
pheno_path = sys.argv[2] #Give the phenotype directory
dataset = sys.argv[3] #nom du jeu de données
param_path = "/home/louison/Documents/FAC/M2/c++_project/detection-of-epistasia-pattern/SMMB-ACO/SMMB-ACO-parameters.txt"
try:
    os.makedirs("SMMB-ACO_results/"+dataset)
    os.makedirs(SMMB-ACO_eval_results)

except FileExistsError:
    print("Directory already exist")

for geno in os.listdir(geno_path):
    geno_dir = "SMMB-ACO_results/"+ dataset + "/" + geno[:-13]
    os.makedirs(geno_dir)
    for i in range(1,21):
        pheno = geno[:-12] + "Phenotype.txt"
        print(geno)
        print(pheno)
        #TODO changer l'ordre des paramètres.
        os.system("./SMMB-ACO/SMMB-ACO "+geno_dir+"/res_"+os.path.basename(geno)+"_"+str(i)+" "+geno_path+"/"+geno+" "+pheno_path+"/"+pheno+ " " +param_path)


#LAUNCH EVAL SCRIPT
#run : ./eval_simu.py input_directory_path output_directory_path n_runs nb_snp len_pattern
#The input directory path must be <Model_results_directory>/<Jeu_donnees_x_directory> with a directory /<fichier_simulé_x> inside containing the n itération for that file.
os.system("./eval/eval_simu.py " "/home/louison/Documents/FAC/M2/c++_project/detection-of-epistasia-pattern/SMMB-ACO_results/"+dataset+ " /home/louison/Documents/FAC/M2/c++_project/detection-of-epistasia-pattern/SMMB-ACO_eval_results 20 100 2")
