/*
 * TypeDef.h
 *
 *  Created on: 19 nov. 2018
 *      Author: courtin
 */

#ifndef TYPEDEF_H_
#define TYPEDEF_H_

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/array.hpp>

typedef boost::numeric::ublas::matrix<int> int_matrix_type;
typedef boost::numeric::ublas::matrix<float> float_matrix_type;
typedef boost::numeric::ublas::matrix<int_matrix_type> matrix_of_int_matrix_type;




#endif /* TYPEDEF_H_ */