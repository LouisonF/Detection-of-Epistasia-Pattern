#! /usr/bin/python3
#-*-coding: utf-8-*-
#pas de shebang pour pouvoir lancer dans l'environnement conda, mettre python devant le ./file_to_launch
# Louison Fresnais M2BB


import sys



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
try:
	prefix = str(args.get("prefix"))
	print("Prefix value: ",prefix)

except TypeError:
	prefix = int(input("enter a disctinctive prefix that will be applied to each result file for this run"))

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



