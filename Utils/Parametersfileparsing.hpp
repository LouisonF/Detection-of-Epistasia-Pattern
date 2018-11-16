/*
 * Parametersfileparsing.hpp
 *
 *  Created on: 16 nov. 2018
 *      Author: louison
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
	unsigned int header_nrows;
	char sep;
	unsigned float alpha;
	unsigned int number_smmbaco_runs;
	unsigned float precision;
	unsigned int number_aco_iter;
	unsigned int number_ants;
	unsigned int number_snp_per_ant;
	unsigned int smallest_subset_size;
	unsigned int max_trials_smmb;
	unsigned int max_trials_learn_mb;
	//evaporation rate parameter
	unsigned int aco_tau_init;
	unsigned int aco_rho;
	unsigned int aco_lambda;
	//Probability distribution parameter
	unsigned int aco_eta;
	unsigned int aco_alpha;
	unsigned int aco_beta;

	unsigned int n_mbs;

	string geno_file;
	string pheno_file;
	//Methods
	Parameters_file_parsing();
	virtual ~Parameters_file_parsing();
	void Parsing();
	void list_parameters() const;
	void import_line(string const& line);
private:
	vector<string> split(string const& s, char delim);
};

#endif /* PARAMETERSFILEPARSING_HPP_ */
