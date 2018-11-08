#! /usr/bin/python3
#-*-coding: utf-8-*-
#pas de shebang pour pouvoir lancer dans l'environnement conda, mettre python devant le ./file_to_launch
# Louison Fresnais M2BB
# François Courtin M2BB

# This script requires numpy (pip install numpy)
import sys
import numpy as np
import re

from parameter_parsing import *


args = vars(argument_parsing())
print(args)


#####################################################################
#   						INPUT	 								#
#####################################################################

###Number of file reading
try:
	number_file = int(args.get("file"))
	print("Number of file value: ",number_file)

except TypeError:
	number_file = int(input("enter the number of file you want to generate in one run"))

###Prefix value reading
prefix = str(args.get('prefix'))
print("Prefix value: ",prefix)
if (prefix == 'None'):

	prefix = str(input("enter a disctinctive prefix that will be applied to each result file for this run"))


###Number of variables reading
try:
	number_variable = int(args.get("variable"))
	print("Number of variable value: ",number_file)

except TypeError:
	number_variable = int(input("enter the number of variables to simulate"))

###Number of cases reading
try:
	number_case = int(args.get("case"))
	print("Number of case value: ",number_case)

except TypeError:
	number_case = int(input("enter the number of cases to simulate"))

###Number of controls reading
try:
	number_control = int(args.get("control"))
	print("Number of control value: ",number_control)

except TypeError:
	number_control = int(input("enter the number of controls to simulate"))

#####################################################################
#   						EXECUTION	 							#
#####################################################################




for i in range(1,number_file+1):
	#file names creation
	geno_file_name = str(prefix+"genotypes"+str(i)+".txt")
	pheno_file_name = str(prefix+"phenotypes"+str(i)+".txt")
	#SNPs IDs matrix creation
	matrix_genotype_ID = []
	for x in range(1,number_variable+1):
		x = str(x)
		ID = "N"+x
		matrix_genotype_ID.append(ID)

	#Convert the ID list into an array
	matrix_genotype_ID = np.asarray(matrix_genotype_ID)
	#Matrix random creation
	matrix_case_geno = np.random.randint(low=0,high=3, size=(number_case, number_variable),dtype="int")
	matrix_control_geno = np.random.randint(low=0,high=3, size=(number_control, number_variable))
	matrix_case_pheno = np.random.randint(low=0,high=2, size=(number_case, 1))
	matrix_control_pheno = np.random.randint(low=0,high=2, size=(number_control, 1))

	#vstack will concatenate by rows case and control array
	matrix_final_geno = np.vstack((matrix_case_geno,matrix_control_geno))
	#save the simulated genotype data file
	np.savetxt(geno_file_name,np.r_[[matrix_genotype_ID], matrix_final_geno], fmt='%s', delimiter=',')

	matrix_final_pheno = np.vstack((matrix_case_pheno,matrix_control_pheno))
	#save the simulated phenotype data file
	np.savetxt(pheno_file_name,np.r_[matrix_final_pheno], fmt='%s', delimiter=',')
