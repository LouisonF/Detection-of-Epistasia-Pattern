/*
 * SmmbACO.hpp
 *
 *  Created on: 9 nov. 2018
 *      Author: Louison Fresnais, François Courtin
 *      Project: SMMB-ACO and Genetic Algorithm for epistasis detection
 *      Under the supervision of Christine Sinoquet(Nantes University)
 *  Modified on: 05 fev 2018
 */

#ifndef SMMBACO_HPP_
#define SMMBACO_HPP_

#include "Parametersfileparsing.hpp"
#include "datainput.hpp"
#include "Miscellaneous.hpp"
#include "Contingency.hpp"
#include "G2_conditional_test_indep.hpp"
#include "G2_test_indep.hpp"

#include <ctime>
#include <cmath>
#include <vector>
#include <libgen.h> //Needed to use basename
#include <list> //Needed for the score map
#include <map>
#include <algorithm>
#include <random>
#include <numeric>
#include <boost/array.hpp>
#include <set>
typedef blas::matrix_column<blas::matrix<int> > blas_column;






using namespace std;

class Smmb_ACO {
public:
	Smmb_ACO(blas::matrix<int> &genos, blas::matrix<int> &phenos , Parameters_file_parsing params);
	virtual ~Smmb_ACO();
	void run_ACO();
	void sum_tau();
	float pheromone_for_snp(float tau_for_snp, float eta_for_snp);
	void compute_distrib_prob();
	void compute_cumulative_distrib_proba();
	unsigned int select_snp_in_distrib_prob(float prob);
	void snp_sampling(vector<unsigned int> &snp_table);
	void learn_mb(list<unsigned int> &mb, vector<unsigned int> &snp_table);
	void forward_phase(list<unsigned int> &mb, vector<unsigned int> &snp_table);
	void backward_phase(list<unsigned int> &mb, vector<unsigned int> &snp_table);
	void evaporation_rate_update(unsigned int snp_index, float g2_score);
	void update_tau();
	void best_mbs(vector<vector<unsigned int>> &mbs);
	void write_results(string output_path);
	map<float, vector<unsigned int>> cumulated_distrib_prob;
private:
	string file_basename;
	string file_path;
	Parameters_file_parsing _params;
	mt19937 rand_seed;
	ofstream _results_handler;
	int number_of_snps;
	blas::matrix<int> &_genotypes;
	blas::matrix<int> &_phenotypes;
	int number_of_indep_test;
	float sum_of_tau;
	vector<float> tau;
	vector<float> eta;
	vector<float> pdf;
	int number_executions;
	ofstream output_file;
	map<unsigned int, vector<double> > scores;
	vector<vector<unsigned int>> mbs;
	map<vector<unsigned int>, vector<double> > results;
	vector<vector<unsigned int>> mbs_count;
};

#endif /* SMMBACO_HPP_ */
