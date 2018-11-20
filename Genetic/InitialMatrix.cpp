/*
 * InitialMatrix.cpp
 *
 *  Created on: 20 nov. 2018
 *      Author: courtin
 */

#include "InitialMatrix.h"

InitialMatrix::InitialMatrix(int_matrix_type Mgeno, int_matrix_type Mpheno) : Mgeno(Mgeno), Mpheno(Mpheno)
{
	// TODO Auto-generated constructor stub

}

InitialMatrix::~InitialMatrix() {
	// TODO Auto-generated destructor stub
}

void InitialMatrix::display_geno(){
	cout << Mgeno << endl;
}

void InitialMatrix::display_pheno(){
	cout << Mpheno << endl;
}
