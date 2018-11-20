from random import *
import math
import numpy as np
import pandas as pd
import itertools


###randrange function but with float
def randrange_float(start, stop, step):
    return randint(0, int((stop - start) / step)) * step + start

#temp
size_epistasia = 3
list_psi_glob = list()
list_random_glob = list()
for logit_iter in range(0,100):

		###Random betas creation
	list_random = list()
	for i in range(0,size_epistasia):
		list_random.append(randrange_float(-1,1,0.1))
	print(list_random)
	list_comb = [[0,1,2]]
	for x in range(1,size_epistasia):
		list_comb.append([0,1,2])
		print(list_comb)
	all_combinations = list(itertools.product(*list_comb))

	print(all_combinations)
	print(len(all_combinations))

	if len(all_combinations) != pow(3,size_epistasia):
		print("Error in combinations")

	list_Psi = list()
	count_healthy = int(0)
	determination_th = (11,12,13,14)
	for i in range(0,pow(3,size_epistasia)):
		temp = all_combinations[i]
		Y = 1
		#multiplicate_Bs = 1
		#multiplicate_Xs = 1
		for i in range(0,len(list_random)):
			Y = Y + list_random[i]*temp[i]
			#We can't iniate multiplicators with 0, that is why at firt iteration, multiplicators take the first value for temp or beta
			if i != 0:
				multiplicate_Bs = multiplicate_Bs*list_random[i] 
				multiplicate_Xs = multiplicate_Xs*temp[i]
			else:
				multiplicate_Bs = list_random[0]
				multiplicate_Xs = temp[0]

		Y =  Y + (multiplicate_Bs*multiplicate_Xs)
		pr = (1/(1+math.exp(-Y)))
		if pr <0.5:
			list_Psi.append(pr)
			count_healthy+=1
	if count_healthy in determination_th:
		list_psi_glob.append(list_Psi)
		list_random_glob.append(list_random)

	#Print list of result( debugging, will print only usefull data in the future)
print("list of the second iteration psi values",list_psi_glob)
print("list of the betas randomly generated in the second iteration",list_random_glob)
#Here we obtain among other informations, a list of Betas that allow us to generate phenotype data accorind to epistasis
#Next step is to do this script as a function and return the list ob betas.