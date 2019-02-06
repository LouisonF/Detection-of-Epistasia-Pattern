/*
 * main.cpp
 *
 *  Created on: 4 déc. 2018
 *      Author: louison
 */



#include <iostream>
#include <stdlib.h>
#include <ctime>
#include <chrono>

#include "SmmbACO.hpp"
#include "common.h"

using namespace std;

int main(int argc, char *argv[])
{

    if(argc < 4)
    {
        cerr << "Missing parameter :\n"
             << "\t./executable <output_path> <path_to_genotypes> <path_to_phenotypes> <param_path>"
             << endl;
        exit(-1);
    }

    // Arguments
    string output = argv[1];
    string genos_file = argv[2];
    string phenos_file = argv[3];
    string param_file = argv[4];
    //string genos_file = "/home/louison/Documents/FAC/M2/c++_project/detection-of-epistasia-pattern/SMMB-ACO/Debug/Naif_1_Genotype.txt";
    //string phenos_file = "/home/louison/Documents/FAC/M2/c++_project/detection-of-epistasia-pattern/SMMB-ACO/Debug/Naif_1_Phenotype.txt";
    //string param_file = "/home/louison/Documents/FAC/M2/c++_project/detection-of-epistasia-pattern/SMMB-ACO/SMMB-ACO-parameters.txt";

//  PARAMETERS
    Parameters_file_parsing params(param_file);
    params.list_parameters();

    params.genos_file = genos_file;      // genotype file name (path)
    params.phenos_file = phenos_file;    // phenotype file name (path)
    int header = params.header_nrows;
    char separator = params.sep;

    //  DATA IMPORTATION
    //Calling Genotype input
    string filename = params.genos_file;
    char sep = params.sep;
    unsigned int header_nrows = params.header_nrows;

    Data_input gen_inst(filename, sep, header_nrows);
    blas::matrix<int> genos = gen_inst.read();

	//Calling Phenotype input
	filename = params.phenos_file;
	Data_input phen_inst(filename, params.sep,params.header_nrows);
	blas::matrix<int> phenos = phen_inst.read();
	//blas_column pheno_col = blas_column(phenos,1);

    cout << endl << "Data imported : " << genos.size1() << " individuals X " << genos.size2() << " SNPs" << endl;

    cout << "pheno size = " << phenos.size1() << endl;
    chrono::high_resolution_clock::time_point t1 = chrono::high_resolution_clock::now();
    Smmb_ACO smmb_ACO(genos, phenos, params);
    smmb_ACO.run_ACO();
    //smmb_ACO.best_mbs(smmb_ACO.mbs);
    smmb_ACO.write_results(output);
    chrono::high_resolution_clock::time_point t2 = chrono::high_resolution_clock::now();
    double duration = chrono::duration_cast<chrono::milliseconds>(t2-t1).count();
    cout << "run time" << duration<<endl;
    cout << "END OF THE PROGRAM"<<endl;

    return 0;
}

