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
#include "SmmbACO.hpp"

#include <iostream>
#include <string>
#include <fstream>
#include <map>
#include <boost/numeric/ublas/matrix.hpp>

using namespace std;

int main()
{
	/*string filename = "test_data.txt";
	unsigned int header_nrows = 2;
	char sep = ',';
	Data_input test(filename, sep, header_nrows);
	test.read();*/
	string file_path = "/home/louison/Documents/FAC/M2/c++_project/detection-of-epistasia-pattern/SMMB-ACO/SMMB-ACO-parameters.txt";
	Parameters_file_parsing test_param(file_path);
	//test_param.Parsing();
	//test_param.list_parameters();

	Smmb_ACO test_smmb(test_param);
	test_smmb.run_ACO();


}

