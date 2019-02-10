/*
 * Parametersfileparsing.hpp
 *
 *  Created on: 16 nov. 2018
 *      Author: Louison Fresnais, François Courtin M2BB
 *      Project: SMMB-ACO and Genetic Algorithm for epistasis detection
 *      Under the supervision of Christine Sinoquet(Nantes University)
 *      Most of this code is from the SMMB-ACO implementation from Clément Niel
 *   Modified on: 19 nov 2018
 */

#ifndef PARAMETERSFILEPARSING_HPP_
#define PARAMETERSFILEPARSING_HPP_

#include <string>
#include <vector>

using namespace std;

class Parameters_file_parsing {
public:
	//Attributes are reachable from any other class
	//Attributes
	string file_path;
	unsigned int header_nrows;
	char sep;
	int len_pop;
	int len_pattern;
	int nb_parents;
	int nb_it;


	//Probability distribution parameter
	float alpha;
	float P_mutation;
	float P_selection;

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
