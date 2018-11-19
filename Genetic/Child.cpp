/*
 * Child.cpp
 *
 *  Created on: 19 nov. 2018
 *      Author: courtin
 */
#include "Child.h"

Child::Child(matrix_of_int_matrix_type MParents_geno, matrix_of_int_matrix_type MParents_pheno, int_matrix_type MParents) : MParents_geno(MParents_geno), MParents_pheno(MParents_pheno), MParents(MParents)
{


}

Child::~Child() {
	// TODO Auto-generated destructor stub
}

void Child::set_children(){

	//Hasard cuts in genotype of parents
	cut1 = rand()%(MParents_geno(0,0).size2()-1)+1;
	cut2 = rand()%(MParents_geno(0,1).size2()-1)+1;
	cut3 = rand()%(MParents_geno(0,2).size2()-1)+1;
	cut4 = rand()%(MParents_geno(0,3).size2()-1)+1;

	//Initialization of 4 children with the good matrix size for the pattern after the cross over
	Mchild1(MParents_geno(0,0).size1(), cut1 + MParents_geno(0,1).size2() - cut2);
	Mchild2(MParents_geno(0,1).size1(), cut2 + MParents_geno(0,0).size2() - cut1);
	Mchild3(MParents_geno(0,2).size1(), cut3 + MParents_geno(0,2).size2() - cut4);
	Mchild4(MParents_geno(0,3).size1(), cut4 + MParents_geno(0,2).size2() - cut3);

}

void Child::start_crossing_over(){
	//Settinf first part of the children genotype
	for (int i = 0; i < int(Mchild1.size1()); i++){
		for (int j = 0; j < cut1; j++){
			Mchild1(i, j) = MParents_geno(0,0)(i, j);
		}
	}

	for (int i = 0; i < int(Mchild2.size1()); i++){
		for (int j = 0; j < cut2; j++){
			Mchild2(i, j) = MParents_geno(0,1)(i, j);
		}
	}

	for (int i = 0; i < int(Mchild3.size1()); i++){
		for (int j = 0; j < cut3; j++){
			Mchild3(i, j) = MParents_geno(0,2)(i, j);
		}
	}

	for (int i = 0; i < int(Mchild4.size1()); i++){
		for (int j = 0; j < cut4; j++){
			Mchild4(i, j) = MParents_geno(0,3)(i, j);
		}
	}
}

