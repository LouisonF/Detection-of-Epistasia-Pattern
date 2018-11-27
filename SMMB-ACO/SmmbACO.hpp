/*
 * SmmbACO.hpp
 *
 *  Created on: 24 nov. 2018
 *      Author: louison
 */

#ifndef SMMBACO_HPP_
#define SMMBACO_HPP_

#include "Parametersfileparsing.hpp"
#include "datainput.hpp"
#include <ctime>
#include <cmath>
#include <vector>
#include <libgen.h> //Needed to use basename




using namespace std;

class Smmb_ACO {
public:
	Smmb_ACO(blas::matrix<int> &genos, blas::matrix<int> &phenos , Parameters_file_parsing params);
	virtual ~Smmb_ACO();
	void run_ACO();
	void sum_tau();
private:
	string file_path;
	Parameters_file_parsing _params;
	int rand_seed;
	ofstream _results_handler;
	int number_of_snps;
	blas::matrix<int> _genotypes;
	blas::matrix<int> _phenotypes;
	int number_of_indep_test;
	float sum_of_tau;
	vector<float> tau;
	vector<float> eta;
	vector<float> pdf;
	int number_executions;
	ofstream output_file;
};

#endif /* SMMBACO_HPP_ */
