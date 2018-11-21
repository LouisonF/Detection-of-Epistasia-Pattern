/*
 * InitialMatrix.h
 *
 *  Created on: 20 nov. 2018
 *      Author: courtin
 */

#ifndef INITIALMATRIX_H_
#define INITIALMATRIX_H_

#include "TypeDef.h"
using namespace std;

class InitialMatrix {
public:
	int_matrix_type Mgeno;
	int_matrix_type Mpheno;
public:
	InitialMatrix(int_matrix_type, int_matrix_type);
	virtual ~InitialMatrix();
	virtual void display_geno();
	virtual void display_pheno();
};

#endif /* INITIALMATRIX_H_ */
