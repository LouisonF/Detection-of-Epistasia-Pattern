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

	double_matrix_type Mpop_geno;
	int_matrix_type MParents;
	int_matrix_type MChildren;
	int len_pattern;
	double P_mutation;
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
