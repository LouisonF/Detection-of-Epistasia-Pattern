/*
 * Child.h
 *
 *  Created on: 19 nov. 2018
 *      Author: Louison Fresnais, François Courtin
 *      Project: SMMB-ACO and Genetic Algorithm for epistasis detection
 *      Under the supervision of Christine Sinoquet(Nantes University)
 */

#ifndef CHILD_H_
#define CHILD_H_
#include "TypeDef.h"
using namespace std;

class Child{
private:

	double_matrix_type Mpop_geno; //Matrix of index of solution of the population in the global matrix
	int_matrix_type MParents; //Matrix of the selected parents
	int_matrix_type MChildren; //Matrix of the children created by crossing over
	int len_pattern;
	double P_mutation; //Probability of mutation
	int nb_snp;

public:
	Child(double_matrix_type, int_matrix_type, int, double, int);
	~Child();
	void set_children();
	void mutation();
	void display_children();
	int_matrix_type get_MChildren();
};

#endif /* CHILD_H_ */
