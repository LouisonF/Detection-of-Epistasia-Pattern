#! /usr/bin/python3
#-*-coding: utf-8-*-
# launch example: python launch_smmb.py ../Simu_naive/Simu_naive_2snp_0.5_geno/ Simu_naive_2snp_0.5
# launch exemple: python launch_smmb.py ../simu_naive_CARLUER_OUEDRAOGO/simu_naive_2snp_0_059_0_25 simu_naive_2snp_0_059_0_25
# Louison Fresnais M2BB
# François Courtin M2BB

import os, os.path
import sys
import time
from subprocess import Popen, PIPE


data_path = sys.argv[1] #Give the data directory
dataset = sys.argv[2] #nom du jeu de données
param_path = sys.argv[3] #path to the parameters file
number_snps = sys.argv[4] # the number of snps in the dataset
size_epistasia = sys.argv[5] # the size of the epistasia pattern

nbr_runs = 2 #change this value in order to modifiy the number of runs
try:
    os.makedirs("SMMB-ACO_results/"+dataset)
    os.makedirs("SMMB-ACO_eval_results")

except FileExistsError:
    print("Directory already exist")
pheno_list = []
pipe_pheno = Popen("ls "+data_path+" | grep -i pheno", shell=True, stdout=PIPE)
for pheno in pipe_pheno.stdout:
    temp_pheno = str(pheno,'utf-8')
    temp_pheno = temp_pheno.strip('\n')
    pheno_list.append(temp_pheno)
j = 0
print(pheno_list)
file = open("time.txt", "a")
start_time = int(round(time.time() * 1000))

pipe_geno = Popen("ls "+data_path+" | grep -i geno", shell=True, stdout=PIPE)
for geno in pipe_geno.stdout:
    geno_temp = str(geno,'utf-8')
    geno_temp = geno_temp.strip('\n')
    geno_temp = geno_temp.strip('.txt')
    print(geno_temp)
    geno_dir = "SMMB-ACO_results/"+dataset+"/"+geno_temp
    try:
        os.makedirs(geno_dir)
    except FileExistsError:
        print("Directory already exist")

    for i in range(1,nbr_runs):
        pheno_str = str(pheno_list[j])
        print(pheno_str)
        exec_start_time = int(round(time.time() * 1000))
        os.system("./SMMB-ACO/smmb_aco "+geno_dir+"/res_"+geno_temp+"_"+str(i)+" "+data_path+"/"+geno_temp+".txt "+data_path+"/"+pheno_str+ " " +param_path)
        exec_end_time = str(float(int(round(time.time() * 1000)) - exec_start_time)/1000)
        file.write(exec_end_time+"s\n")

    j+=1
end_time = str(float((int(round(time.time() * 1000)) - start_time))/1000/60)
print(end_time)
file.close()

#LAUNCH EVAL SCRIPT
#run : ./eval_simu.py input_directory_path output_directory_path n_runs nb_snp len_pattern
#The input directory path must be <Model_results_directory>/<Jeu_donnees_x_directory> with a directory /<fichier_simulé_x> inside containing the n itération for that file.
os.system("./eval/eval_simu.py " "/comptes/E146938Q/c++_project/rendu_smmb_aco/SMMB-ACO_results/"+dataset+ " /comptes/E146938Q/c++_project/rendu_smmb_aco/SMMB-ACO_eval_results "+str(nbr_runs)+" "+str(number_snps)+" "+str(size_epistasia))
