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
	int_matrix_type Mgeno;
	int_matrix_type Mpheno;
	int len_pop;
	int len_patern;
	int nb_sol;
	int nb_indiv;
	int_matrix_type select_sol;
	matrix_of_int_matrix_type Mpop_geno;
	matrix_of_int_matrix_type Mpop_pheno;
public:
	Population(int_matrix_type, int_matrix_type, int, int);
	virtual ~Population();
	void sol_selection();
	void init_pop_geno();
	void init_pop_pheno();
	void display_geno_sol();
	void display_pheno_sol();
};

#endif /* POPULATION_H_ */
