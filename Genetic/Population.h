/*
 * Population.h
 *
 *  Created on: 19 nov. 2018
 *      Author: courtin
 */

#ifndef POPULATION_H_
#define POPULATION_H_

#include "TypeDef.h"
#include "InitialMatrix.h"

class Population : public InitialMatrix {
protected:
	int nb_sol;
	int len_pattern;
	int_matrix_type Mpop_geno;
public:
	Population(int_matrix_type, int_matrix_type, int, int);
	~Population();
	void sol_selection();
	void init_pop_geno();
	void display_geno_sol();
	int_matrix_type get_Mpop_geno();
};

#endif /* POPULATION_H_ */
