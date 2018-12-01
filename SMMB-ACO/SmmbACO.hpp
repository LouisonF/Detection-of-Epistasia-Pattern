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
#include <list> //Needed for the score map
#include <map>
#include <algorithm>
#include <random>;





using namespace std;

class Smmb_ACO {
public:
	Smmb_ACO(blas::matrix<int> &genos, blas::matrix<int> &phenos , Parameters_file_parsing params);
	virtual ~Smmb_ACO();
	void run_ACO();
	void sum_tau();
	float pheromone_for_snp(float tau_for_snp, float eta_for_snp);
	void compute_distrib_prob();
	void compute_cumlative_dristrib_proba();
	unsigned int select_snp_in_distrib_prob(float prob);
	void snp_sampling(vector<unsigned int> &snp_table);
	void forward_phase(vector<unsigned int> &mb, vector<unsigned int> &snp_table);
	//void backward_pahse(vector<unsigned int> &mb, vector<unsigned int> &snp_table)

private:
	string file_path;
	Parameters_file_parsing _params;
	mt19937 rand_seed;
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
	map<unsigned, list<float> > scores;
	map<float, vector<unsigned>> cumulated_distrib_prob;
};

#endif /* SMMBACO_HPP_ */
