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

using namespace std;

class Smmb_ACO {
public:
	Smmb_ACO();
	virtual ~Smmb_ACO();
	void run_aco();
private:
	string file_path;
};

#endif /* SMMBACO_HPP_ */
