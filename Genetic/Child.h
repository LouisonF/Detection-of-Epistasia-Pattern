/*
 * Child.h
 *
 *  Created on: 19 nov. 2018
 *      Author: courtin
 */

#ifndef CHILD_H_
#define CHILD_H_
#include "Population.h"

class Child {
public:

	matrix_of_int_matrix_type MParents_geno;
	matrix_of_int_matrix_type MParents_pheno;

	int_matrix_type MParents;

	matrix_of_int_matrix_type MChild_geno;
	matrix_of_int_matrix_type MChild_pheno;

	int cut1;
	int cut2;
	int cut3;
	int cut4;

	int_matrix_type Mchild1;
	int_matrix_type Mchild2;
	int_matrix_type Mchild3;
	int_matrix_type Mchild4;

	Child(matrix_of_int_matrix_type, matrix_of_int_matrix_type, int_matrix_type);
	~Child();
	void set_children();
	void start_crossing_over();
	void end_crossing_over();
};

#endif /* CHILD_H_ */
