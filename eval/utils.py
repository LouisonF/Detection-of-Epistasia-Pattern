#! /usr/bin/python3
#-*-coding: utf-8-*-
# In order to launch the script in a conda environment, use this call: " Python ./file_to_launch"
# Louison Fresnais M2BB
# François Courtin M2BB

def write_res(output_path, value):

    path = output_path[0:(len(output_path)-4)]
    file = open(output_path+"_result.txt", "a")
    file.write(value+"\n")
    file.close()


#************************************************************


def create_res_files(input, output, caus_snp):

    caus_snp_array = caus_snp.split(',')

    TP = 0
    FP = 0
    FN = 0


    file_list = glob.glob(input + "/*.txt")

    for files in file_list:
        try:
            output_path = output + "/" + os.path.basename(files)
        except AttributeError:
            print('Error in the input and output names')

        file  = open(files, "r")
        #print(file.read())

        for line in file:
            pres_TP = False
            inc_TP = False
            inc_FP = False
            inc_FN = False
            try:
                pattern = re.search('{(.+?)}', line).group(1)
            except AttributeError:
                pattern = ''
            if pattern == '':
                write_res(output_path, "FN")
            else:
                snp_array = pattern.split(",")
                if (snp_array == caus_snp_array) or (snp_array[::-1] == caus_snp_array):
                    write_res(output_path, "TP")
                elif len(snp_array) < len(caus_snp_array):
                    write_res(output_path, "FP")
                else:
                    for i in range(len(snp_array)):
                        for j in range(len(caus_snp_array)):
                            if snp_array[i] == caus_snp_array[j]:
                                pres_TP = True
                if (pres_TP):
                    write_res(output_path, "TP")
                else:
                    write_res(output_path, "FP")
