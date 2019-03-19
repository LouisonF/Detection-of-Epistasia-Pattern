/*
 * Parent.h
 *
 *  Created on: 19 nov. 2018
 *      Author: Louison Fresnais, François Courtin
 *      Project: SMMB-ACO and Genetic Algorithm for epistasis detection
 *      Under the supervision of Christine Sinoquet(Nantes University)
 */

#ifndef PARENT_H_
#define PARENT_H_
#include "TypeDef.h"

class Parent{

protected:
	int nb_sol;
	int nb_parents;
	int len_pattern;
	double median;
	double P_selection;
	double_matrix_type Mpop_geno;
	int_matrix_type Mparents;

public :
	Parent(int, int, int, double, double, double_matrix_type);
	~Parent();
	void parents_selection();
	int_matrix_type get_MParents();
	void display_parents();
};

#endif /* PARENT_H_ */
