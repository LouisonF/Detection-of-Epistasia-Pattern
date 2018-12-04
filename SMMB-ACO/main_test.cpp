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
#include "Parametersfileparsing.hpp"*
#include "Miscellaneous.hpp"
#include "SmmbACO.hpp"

#include <iostream>
#include <string>
#include <fstream>
#include <map>
#include <boost/numeric/ublas/matrix.hpp>

using namespace std;

int douze()
{
	/*//Calling parameters file parsing
	string file_path = "/home/louison/Documents/FAC/M2/c++_project/detection-of-epistasia-pattern/SMMB-ACO/SMMB-ACO-parameters.txt";
	Parameters_file_parsing test_param(file_path);
	test_param.Parsing();

	string filename;
	char sep;
	unsigned int header_nrows;
	//Calling Genotype input
	 filename = test_param.genos_file;
	 sep = test_param.sep;
	 header_nrows = test_param.header_nrows;


	Data_input gen_inst(filename, sep, header_nrows);
	blas::matrix<int> genos = gen_inst.read();

	//Calling Phenotype input
	filename = test_param.phenos_file;
	Data_input phen_inst(filename, test_param.sep,test_param.header_nrows);
	blas::matrix<int> phenos = phen_inst.read();


	Smmb_ACO test_smmb(genos, phenos, test_param);
	test_smmb.run_ACO();
*/
	unsigned int size = 5;
	vector<unsigned int> snp_sorted;
	snp_sorted.push_back(1);
	snp_sorted.push_back(2);
	snp_sorted.push_back(3);
	snp_sorted.push_back(4);
	snp_sorted.push_back(5);
	vector<vector<unsigned int>> all_snps_combinations;
	//Miscellaneous::powerset(snp_sorted,all_index_combinations);
	//Miscellaneous::print_powerset(all_index_combinations);
	//
	Miscellaneous::combinator(snp_sorted,all_snps_combinations,size);
	Miscellaneous::link_comb_to_snp(snp_sorted,all_snps_combinations);

	Miscellaneous::print_human_readable_combinations(all_snps_combinations);
	cout << all_snps_combinations.size();
	return 0;


}

