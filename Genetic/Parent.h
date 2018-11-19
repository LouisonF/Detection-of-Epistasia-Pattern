/*
 * Parent.h
 *
 *  Created on: 19 nov. 2018
 *      Author: courtin
 */

#ifndef PARENT_H_
#define PARENT_H_
#include "Population.h"

class Parent {
public:

	matrix_of_int_matrix_type Mpop_geno;
	matrix_of_int_matrix_type Mpop_pheno;
	int_matrix_type Mparents;
	matrix_of_int_matrix_type MParents_geno;
	matrix_of_int_matrix_type MParents_pheno;


	Parent(matrix_of_int_matrix_type, matrix_of_int_matrix_type);
	~Parent();
};

#endif /* PARENT_H_ */
