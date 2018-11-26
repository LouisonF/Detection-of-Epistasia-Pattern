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
#include <libgen.h> //Needed to use basename




using namespace std;

class Smmb_ACO {
public:
	Smmb_ACO(Parameters_file_parsing params);
	virtual ~Smmb_ACO();
	void run_ACO();
private:
	string file_path;
	Parameters_file_parsing _params;
	int rand_seed;
	ofstream _results_handler;
};

#endif /* SMMBACO_HPP_ */
