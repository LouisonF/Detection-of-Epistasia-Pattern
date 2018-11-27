/*
 * Parent.h
 *
 *  Created on: 19 nov. 2018
 *      Author: courtin
 */

#ifndef PARENT_H_
#define PARENT_H_
#include "TypeDef.h"

class Parent{

protected:
	int nb_sol;
	int nb_parents;
	int len_pattern;
	float median;
	float P_selection;
	float_matrix_type Mpop_geno;
	int_matrix_type Mparents;

public :
	Parent(int, int, int, float, float, float_matrix_type);
	~Parent();
	void parents_selection();
	int_matrix_type get_MParents();
	void display_parents();
};

#endif /* PARENT_H_ */
