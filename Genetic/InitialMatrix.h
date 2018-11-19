/*
 * InitialMatrix.h
 *
 *  Created on: 19 nov. 2018
 *      Author: courtin
 */

#ifndef INITIALMATRIX_H_
#define INITIALMATRIX_H_

#include "TypeDef.h"

class InitialMatrix {
private:
	int_matrix_type Mgeno;
	int_matrix_type Mpheno;
public:
	InitialMatrix();
	virtual ~InitialMatrix();

};

#endif /* INITIALMATRIX_H_ */
