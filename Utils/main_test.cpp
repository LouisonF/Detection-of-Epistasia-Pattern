/*
 * main_test.cpp
 *
 *  Created on: 13 nov. 2018
 *      Author: Louison Fresnais, François Courtin M2BB
 *      Project: SMMB-ACO and Genetic Algorithm for epistasis detection
 *      Under the supervision of Christine Sinoquet(Nantes University)
 *  Modified on: 19 nov 2018
 */
#include "datainput.hpp"
#include "Parametersfileparsing.hpp"

#include <iostream>
#include <string>
#include <fstream>
#include <map>
#include <boost/numeric/ublas/matrix.hpp>

using namespace std;

int main()
{
	string filename = "ouigenotypes1.txt";
	unsigned int header_nrows = 1;
	char sep = ',';
	Data_input test(filename, sep, header_nrows);
	blas::matrix<int> matrice = test.read();
	vector<string>  list = test.get_snps();

    /*for(unsigned i=0;i<matrice.size1();++i)
    {
        cout<<"| ";
        for (unsigned j=0;j<matrice.size2();++j)
        {
            cout<<matrice(i,j)<<" | ";
        }
        cout<<"|"<<endl;
    }*/
	string file_path = "/home/louison/Documents/FAC/M2/c++_project/detection-of-epistasia-pattern/SMMB-ACO/SMMB-ACO-parameters.txt";
	Parameters_file_parsing test_param(file_path);
	test_param.Parsing();
	test_param.list_parameters();

}

