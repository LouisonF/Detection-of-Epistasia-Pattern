/*
 * Parent.h
 *
 *  Created on: 19 nov. 2018
 *      Author: courtin
 */

#ifndef PARENT_H_
#define PARENT_H_
#include "Population.h"

class Parent : public Population{
protected:
	int_matrix_type Mparents;
	matrix_of_int_matrix_type Mparents_geno;
	matrix_of_int_matrix_type Mparents_pheno;

public:
	Parent(matrix_of_int_matrix_type);
	virtual ~Parent();
};

#endif /* PARENT_H_ */
