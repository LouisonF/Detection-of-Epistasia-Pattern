/*
 * Parametersfileparsing.hpp
 *
 *  Created on: 9 nov. 2018
 *      Author: Louison Fresnais, François Courtin M2BB
 *      Project: SMMB-ACO and Genetic Algorithm for epistasis detection
 *      Under the supervision of Christine Sinoquet(Nantes University)
 *      Most of this code is from the SMMB-ACO implementation from Clément Niel
 *   Modified on: 05 fev 2019
 */

#ifndef PARAMETERSFILEPARSING_HPP_
#define PARAMETERSFILEPARSING_HPP_

#include <string>
#include <vector>

using namespace std;

class Parameters_file_parsing {
public:
	//Attributes are reachable from any other class
	//parameters explanation are in the parameter file of SMMB-ACO
	//Attributes
	string file_path;
	unsigned int header_nrows;
	char sep;
	float alpha;
	unsigned int number_smmbaco_runs;
	float precision;
	unsigned int aco_set_size;
	unsigned int number_aco_iter;
	unsigned int number_ants;
	unsigned int number_snp_per_ant;
	unsigned int smallest_subset_size;
	unsigned int size;
	unsigned int max_trials_smmb;
	unsigned int max_trials_learn_mb;
	//evaporation rate parameter
	float aco_tau_init;
	float aco_rho;
	float aco_lambda;
	//Probability distribution parameter
	float aco_eta;
	float aco_alpha;
	float aco_beta;

	unsigned int n_mbs;

	string genos_file;
	string phenos_file;
	//Methods
	Parameters_file_parsing(const string);
	Parameters_file_parsing();
	virtual ~Parameters_file_parsing();
	void Parsing();
	void list_parameters();
	void import_line(string const line);
private:
	vector<string> split(string const s, char delim);
};

#endif /* PARAMETERSFILEPARSING_HPP_ */
