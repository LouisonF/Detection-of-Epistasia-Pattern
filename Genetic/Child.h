/*
 * Child.h
 *
 *  Created on: 19 nov. 2018
 *      Author: courtin
 */

#ifndef CHILD_H_
#define CHILD_H_
#include "TypeDef.h"
using namespace std;

class Child{
private:

	float_matrix_type Mpop_geno;
	int_matrix_type MParents;
	int_matrix_type MChildren;
	int len_pattern;
	float P_mutation;
	int nb_snp;

public:
	Child(float_matrix_type, int_matrix_type, int, float, int);
	~Child();
	void set_children();
	void mutation();
	void display_children();
	int_matrix_type get_MChildren();
};

#endif /* CHILD_H_ */
