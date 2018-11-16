from random import *
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import sklearn
from sklearn.linear_model import LogisticRegression
#temp
size_epistasia = 2
list_psi_glob = list()
list_random_glob = list()
for x in range(1,10):

	###Random betas creation
	list_random = list()
	for i in range(0,size_epistasia):
		list_random.append(random())
	print(list_random)

	###Epistasia matrix
	matrix_epistasia = np.array([[0,0]])
	for x in range(0,3):
		for y in range(0,3):
			temp = np.array([[x,y]])
			matrix_epistasia = np.concatenate((matrix_epistasia,temp),axis=0)
	matrix_epistasia = np.delete(matrix_epistasia,0,0)

	j=0
	list_Psi = list()
	for i in range(0,9):
	 	temp = (matrix_epistasia[i,])
	 	temp = temp.tolist()
	 	print(temp)
	 	Y = list_random[0]*temp[0]+list_random[1]*temp[1]
	 	print("Psi value = ",Y)
	 	list_Psi.append(Y)

	 	j+=1
	list_psi_glob.append(list_Psi)
	list_random_glob.append(list_random)
	 #TODO: Ici on ne fait pas réellement une regression logistique ... comment va-t-on prédire les phénotypes ?



print("list of the second iteration psi values",list_psi_glob[2])
print("list of the betas randomly generated in the second iteration",list_random_glob[2])
