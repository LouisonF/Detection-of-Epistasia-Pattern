#! /usr/bin/python3
#-*-coding: utf-8-*-

#Launch exemple :
# ./launch_genet.py /home/courtin/Documents/M2/ProjetC/smu_naive_100/Simu_naive_2snp_0.5 naive_2snp_05 /home/courtin/Documents/M2/ProjetC/detection-of-epistasia-pattern/Genetic/PARAMETERS_GENETIC.txt 1 100 2 20
# Louison Fresnais M2BB
# François Courtin M2BB

import os, os.path
import sys
import time
from subprocess import Popen, PIPE


data_path = sys.argv[1]
dataset = sys.argv[2]
param_path = sys.argv[3]
n_runs = sys.argv[4]
nb_snp = sys.argv[5]
len_pattern = sys.argv[6]

try:
    os.makedirs("Genetic/Genetic_results/"+dataset)
    os.makedirs("Genetic/Genetic_eval_results")

except FileExistsError:
    print("Directory already exists")

pheno_list = []
pipe_pheno = Popen("ls "+data_path+" | grep -i pheno", shell=True, stdout=PIPE)
for pheno in pipe_pheno.stdout:
    temp_pheno = str(pheno,'utf-8')
    temp_pheno = temp_pheno.strip('\n')
    pheno_list.append(temp_pheno)
j = 0
#print(pheno_list)
file = open("time_"+dataset+".txt", "a")

pipe_geno = Popen("ls "+data_path+" | grep -i geno", shell=True, stdout=PIPE)
for geno in pipe_geno.stdout:
    geno_temp = str(geno,'utf-8')
    geno_temp = geno_temp.strip('\n')
    geno_temp = geno_temp.strip('.txt')
    #print(geno_temp)
    geno_dir = "Genetic/Genetic_results/"+dataset+"/"+geno_temp
    try:
        os.makedirs(geno_dir)
    except FileExistsError:
        print("Directory already exists")

    print("Run Genetic on " + geno_temp)
    for i in range(1,int(n_runs)+1):
        pheno_str = str(pheno_list[j])
        #print(pheno_str)
        exec_start_time = int(round(time.time() * 1000))
        os.system("./Genetic/Release/Genetic "+geno_dir+"/res_"+geno_temp+"_"+str(i)+" "+data_path+"/"+geno_temp+".txt "+data_path+"/"+pheno_str +" "+param_path)
        exec_end_time = str(float(int(round(time.time() * 1000)) - exec_start_time)/1000)
        file.write(exec_end_time+"s\n")

    j+=1



#LAUNCH EVAL SCRIPT
#run : ./eval_simu.py input_directory_path output_directory_path n_runs nb_snp len_pattern
#The input directory path must be <Model_results_directory>/<Jeu_donnees_x_directory> with a directory /<fichier_simulé_x> inside containing the n itération for that file.
os.system("./eval/eval_simu.py " "Genetic/Genetic_results/"+dataset+ " Genetic/Genetic_eval_results "+n_runs+" "+nb_snp+" "+len_pattern)
