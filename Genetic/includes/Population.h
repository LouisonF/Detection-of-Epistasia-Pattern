/*
 * Population.h
 *
 *  Created on: 19 nov. 2018
 *      Author: Louison Fresnais, François Courtin
 *      Project: SMMB-ACO and Genetic Algorithm for epistasis detection
 *      Under the supervision of Christine Sinoquet(Nantes University)
 */

#ifndef POPULATION_H_
#define POPULATION_H_

#include "TypeDef.h"


class Population{
protected:
	int_matrix_type Mgeno;
	int_matrix_type Mpheno;
	int nb_sol;
	int len_pattern;
	double_matrix_type Mpop_geno;
public:
	Population(int_matrix_type, int_matrix_type, int &len_pop, int);
	~Population();
	long double fact(long int N);
	void check_pop_len();
	void sol_selection();
	void init_pop_geno();
	void display_geno_sol();
	double_matrix_type get_Mpop_geno();
	void set_Mpop_geno(int, int, double);
};

#endif /* POPULATION_H_ */
