/*
 * TypeDef.h
 *
 *  Created on: 19 nov. 2018
 *      Author: Louison Fresnais, François Courtin
 *      Project: SMMB-ACO and Genetic Algorithm for epistasis detection
 *      Under the supervision of Christine Sinoquet(Nantes University)
 */

#ifndef TYPEDEF_H_
#define TYPEDEF_H_

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/array.hpp>

typedef boost::numeric::ublas::matrix<int> int_matrix_type;
typedef boost::numeric::ublas::matrix<double> double_matrix_type;
typedef boost::numeric::ublas::matrix<int_matrix_type> matrix_of_int_matrix_type;




#endif /* TYPEDEF_H_ */
