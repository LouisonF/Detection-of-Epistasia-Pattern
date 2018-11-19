/*
 * Population.h
 *
 *  Created on: 19 nov. 2018
 *      Author: courtin
 */

#ifndef POPULATION_H_
#define POPULATION_H_

#include "TypeDef.h"

class Population {
protected:
	int nb_sol;
	int nb_indiv;
	int_matrix_type select_sol;
	matrix_of_int_matrix_type Mpop_geno;
	int_matrix_type Mpop_pheno;
public:
	Population(matrix_of_int_matrix_type, int_matrix_type, int, int);
	virtual ~Population();
	void display_geno_sol();
	void display_pheno_sol();
};

#endif /* POPULATION_H_ */
