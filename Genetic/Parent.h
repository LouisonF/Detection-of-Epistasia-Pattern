/*
 * Parent.h
 *
 *  Created on: 19 nov. 2018
 *      Author: courtin
 */

#ifndef PARENT_H_
#define PARENT_H_
#include "Population.h"

class Parent : public Population {

protected:
	int nb_parents;
	int_matrix_type Mparents;

public :
	Parent(int_matrix_type, int_matrix_type, int, int, int);
	~Parent();
	void parents_selection();
	void display_parents();
};

#endif /* PARENT_H_ */
