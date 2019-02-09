/*
 * G2_conditional_test_indep.cpp
 *
 *  Created on: 9 nov. 2018
 *      Author: Louison Fresnais, François Courtin
 *      Project: SMMB-ACO and Genetic Algorithm for epistasis detection
 *      Under the supervision of Christine Sinoquet(Nantes University)
 *  Modified on: 09 fev 2018
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
	//Here we check if there is the required number of arguments. If not program will exit with an error message
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

    /*string output = "/home/louison/Documents/FAC/M2/c++_project/detection-of-epistasia-pattern/SMMB-ACO_results/Naif_1";
    string genos_file = "/home/louison/Documents/FAC/M2/c++_project/detection-of-epistasia-pattern/SMMB-ACO/Debug/Naif_1_Genotype.txt";
    string phenos_file = "/home/louison/Documents/FAC/M2/c++_project/detection-of-epistasia-pattern/SMMB-ACO/Debug/Naif_1_Phenotype.txt";
    string param_file = "/home/louison/Documents/FAC/M2/c++_project/detection-of-epistasia-pattern/SMMB-ACO/SMMB-ACO-parameters_2snp.txt";*/


    //Parameters file parsing and parameters variable initialisation
    Parameters_file_parsing params(param_file);
    params.list_parameters();

    params.genos_file = genos_file;      // genotype file name (path)
    params.phenos_file = phenos_file;    // phenotype file name (path)

    //  DATA IMPORTATION

    //Calling Genotype input
    string filename = params.genos_file;
    char sep = params.sep;
    unsigned int header_nrows = params.header_nrows;
    //Instantiation of data_input object
    Data_input gen_inst(filename, sep, header_nrows);
    //Genotypes matrix inialisation
    blas::matrix<int> genos = gen_inst.read();

	//Calling Phenotype input
	filename = params.phenos_file;
    //Instantiation of data_input object
	Data_input phen_inst(filename, params.sep,params.header_nrows);
    //Phenotypes matrix inialisation
	blas::matrix<int> phenos = phen_inst.read();
	//Output summary of imported datas
    cout << endl << "Data imported : " << genos.size1() << " individuals X " << genos.size2() << " SNPs" << endl;
    //Time handling (just for information if ran without launcher).
    chrono::high_resolution_clock::time_point t1 = chrono::high_resolution_clock::now();
    //Instanciation of a Smmb_ACO object
    Smmb_ACO smmb_ACO(genos, phenos, params);
    //Call the run method that will compute Smmb_ACO on the input dataset
    smmb_ACO.run_ACO();
    //Call the method that will write the results in the proper file. Check path in the launcher and in the method itself for further information.
    smmb_ACO.write_results(output);
    //Compare time at the beginning of the program versus now. Give the running time of the Smmb_ACO iteration.
    chrono::high_resolution_clock::time_point t2 = chrono::high_resolution_clock::now();
    double duration = chrono::duration_cast<chrono::seconds>(t2-t1).count();
    cout << "run time for this iteration" << duration << " in seconds "<<endl;
    cout << "END OF THE PROGRAM"<<endl;

    return 0;
}

