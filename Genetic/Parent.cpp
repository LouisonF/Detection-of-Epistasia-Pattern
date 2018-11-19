/*
 * Parent.cpp
 *
 *  Created on: 19 nov. 2018
 *      Author: courtin
 */

#include "Parent.h"
#include "Population.h"

Parent::Parent(matrix_of_int_matrix_type Mpop_geno, matrix_of_int_matrix_type Mpop_pheno) : Mpop_geno(Mpop_geno), Mpop_pheno(Mpop_pheno)
{

	Mparents(1, 4);
	//Hasard selection of 2 pairs of parents
	Mparents(0, 0) = rand()%(Mpop_geno.size1());
	Mparents(0, 1) = rand()%(Mpop_geno.size1());
	Mparents(0, 2) = rand()%(Mpop_geno.size1());
	Mparents(0, 3) = rand()%(Mpop_geno.size1());


	//Control if parents are not the same
	while (Mparents(0, 1) == Mparents(0, 0)){
		Mparents(0, 1) = rand()%(Mpop_geno.size1());
	}

	while (Mparents(0, 2) == Mparents(0, 0) or Mparents(0, 2) == Mparents(0, 1)){
		Mparents(0, 2) = rand()%(Mpop_geno.size1());
	}

	while (Mparents(0, 3) == Mparents(0, 0) or Mparents(0, 3) == Mparents(0, 1) or Mparents(0, 3) == Mparents(0, 2)){
		Mparents(0, 3) = rand()%(Mpop_geno.size1());
	}

	MParents_geno(1, 4);
    for (int i = 0; i < int(MParents_geno.size2()); i++){
    	MParents_geno(0, i) = Mpop_geno(Mparents(0, i), 0);
    }

   MParents_pheno(1, 4);
    for (int i = 0; i < int(MParents_pheno.size2()); i++){
    	MParents_pheno(0, i) = Mpop_pheno(Mparents(0, i), 0);
    }
}

Parent::~Parent() {
	// TODO Auto-generated destructor stub
}

