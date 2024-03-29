from random import *
import math
import itertools


###randrange function but with float

def randrange_float(start, stop, step):
    return randint(0, int((stop - start) / step)) * step + start
###The epistasis_combination function take in parameter the size of epistasis pattern and the list of combinations.
#It will return all possible combinations for a pattern of given size

def epistasis_combination(size_epistasia,list_comb):
    for x in range(1,size_epistasia):
        list_comb.append([0,1,2])
    all_combinations = list(itertools.product(*list_comb))
    return all_combinations
###The random_betas_list function produce a list of randomly created betas between
#  0 and 1 according to the size of the pattern (that will determine the number
#of betas in the logit model)

def random_betas_list(size_epistasia):
    list_random = list()
    for i in range(0,size_epistasia):
        list_random.append(randrange_float(0,1,0.1))
    return list_random
###The determination_th function takes all combinations and an allowed
#error_percentage and return the threshold value that will be used to determine
#if a logit model is good enough to produce causal snps.

def determination_th(all_combinations,error_percentage):
     determination_th = ()
     prcentage = int(len(all_combinations)*error_percentage)
     if prcentage <1:
         prcentage = 1
     for x in range(-prcentage,prcentage):
          determination_th = determination_th + (int(len(all_combinations)/2+x),)
     return determination_th

### The compute_logit function will compute a logit model according to a combination
# and Betas values in the list of betas.
def compute_logit(list_random,combination):
    Y = -1
    for i in range(0,len(list_random)):
      Y = Y + list_random[i]*combination[i]
         #We can't iniate multiplicators with 0, that is why at firt iteration, multiplicators take the first value for temp or beta
      if i != 0:
          multiplicate_Bs = multiplicate_Bs*list_random[i]
          multiplicate_Xs = multiplicate_Xs*combination[i]
      else:
          multiplicate_Bs = list_random[0]
          multiplicate_Xs = combination[0]

    Y =  Y + (multiplicate_Bs*multiplicate_Xs)
    pr = (1/(1+math.exp(-Y)))
    return pr

### THe fit_relevant_logit function is going to call the compute_logit function
#and batas values for a model that allow to produce causal snps for a specified
#pattern size.

def fit_relevant_logit(size_epistasia,all_combinations,determination_th,max_iter):
     list_psi_glob = list()
     list_random_glob = list()
     for logit_iter in range(0,max_iter):
          list_random = random_betas_list(size_epistasia)
          if len(all_combinations) != pow(3,size_epistasia):
              print("Error in combinations")

          list_Psi = list()
          count_healthy = int(0)
          for i in range(0,pow(3,size_epistasia)):
              temp = all_combinations[i]
              pr = compute_logit(list_random,temp)
              if pr <0.5:
                  list_Psi.append(pr)
                  count_healthy+=1
          if count_healthy in determination_th:
               list_psi_glob.append(list_Psi)
               list_random_glob.append(list_random)
     return list_random_glob
#Here we obtain among other informations, a list of Betas that allow us to generate phenotype data according to epistasis.
