#! /usr/bin/python3
#-*-coding: utf-8-*-
# In order to launch the script in a conda environment, use this call: " Python ./file_to_launch"
# Louison Fresnais M2BB
# François Courtin M2BB
import re

#write value in the output file
def write_res(output, value, res_type):

    #for TP FP or FN
    if res_type == "result":
        file = open(output+"_results.txt", "a")
    #for the f_measure
    elif res_type == "f_measure":
        file = open(output+"/f_measure.txt", "a")
    #for the power
    elif res_type == "power":
        file = open(output+"/power.txt", "a")
    file.write(value+"\n")
    file.close()
    return(file.name)


#************************************************************#

#Set TP of FP or FN for one result file
def set_res(file, caus_snp):

    caus_snp_array = caus_snp.split(',')
    caus_snp_array = caus_snp_array[::-1] #Reverse the array to match the right sort

    res_file  = open(file, "r")
    for line in res_file:
        #print(line)
        if line[0] != "{":
            continue;
        try:
            pattern = re.search('{(.+?)}', line).group(1)
        except AttributeError:
            #print("c'est l'exception12")
            pattern = ''
        if pattern == '':
            return("FN")
        else:
            snp_array = pattern.split(",")
            if (caus_snp in pattern) or (caus_snp in ','.join(snp_array[::-1])):
                return("TP")
            elif len(snp_array) < len(caus_snp_array):
                return("FN")
            # else:
            #     for i in range(len(snp_array)):
            #         for j in range(len(caus_snp_array)):
            #             if snp_array[i] == caus_snp_array[j]:
            #                 return("TP")
        return("FP")
    res_file.close()

#************************************************************#

#calcul the f-measure
def run_f_measure(input):

    TP = 0
    FP = 0
    FN = 0

    file = open(input,'r')

    for line in file:
        if "TP" in line:
            TP += 1
        elif "FP" in line:
            FP += 1
        elif "FN" in line:
            FN += 1

    if (TP == 0 and FP == 0 and FN != 0):
        FP = 1
    if (TP == 0 and FN == 0 and FP != 0):
        FN = 1

    recall = TP/(TP+FN)
    precision = TP/(TP+FP)

    if TP == 0:
        f_measure = 0
    else:
        f_measure = 2/((1/recall)+(1/precision))

    #print("F-measure : ",f_measure)
    file.close()
    return(f_measure)



#************************************************************#

#calcul the f-measure
def run_power(input,n_runs):

    TP = 0

    file = open(input,'r')

    for line in file:
        if "TP" in line:
            TP += 1

    power = TP/n_runs

    #print("Power : ",power)
    file.close()
    return(power)
