/*
 * Parent.cpp
 *
 *  Created on: 19 nov. 2018
 *      Author: courtin
 */

#include "Parent.h"

Parent::Parent(matrix_of_int_matrix_type Mpop) {

	Mparents(1, 4);
	//Hasard selection of 2 pairs of parents
	Mparents(0, 0) = rand()%(Mpop.size1());
	Mparents(0, 1) = rand()%(Mpop.size1());
	Mparents(0, 2) = rand()%(Mpop.size1());
	Mparents(0, 3) = rand()%(Mpop.size1());


	//Control if parents are not the same
	while (Mparents(0, 1) == Mparents(0, 0)){
		Mparents(0, 1) = rand()%(Mpop.size1());
	}

	while (Mparents(0, 2) == Mparents(0, 0) or Mparents(0, 2) == Mparents(0, 1)){
		Mparents(0, 2) = rand()%(Mpop.size1());
	}

	while (Mparents(0, 3) == Mparents(0, 0) or Mparents(0, 3) == Mparents(0, 1) or Mparents(0, 3) == Mparents(0, 2)){
		Mparents(0, 3) = rand()%(Mpop.size1());
	}
}

Parent::~Parent() {
	// TODO Auto-generated destructor stub
}

