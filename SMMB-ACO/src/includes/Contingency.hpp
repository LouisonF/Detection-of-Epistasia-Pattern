/*
 * datainput.h
 *
 *  Created on: 9 nov. 2018
 *      Author: Louison Fresnais, François Courtin
 *      Project: SMMB-ACO and Genetic Algorithm for epistasis detection
 *      Under the supervision of Christine Sinoquet(Nantes University)
 *      Imported from Clément Niel's code
 *  Modified on: 05 fev 2018
 */
#ifndef CONTINGENCY_HPP
#define CONTINGENCY_HPP

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>

#include "common.h"

namespace blas=boost::numeric::ublas;


class Contingency : public blas_dmatrix
{
public:
    Contingency();
    Contingency(int a, int b);
    Contingency(blas_matrix const& m);
    Contingency(Contingency const& m);
    Contingency(blas_column const& phenos,
                blas_matrix const& genos,
                unsigned index);
    Contingency(blas_column const& phenos, blas_column const& genos);
    Contingency(blas_vector const& phenos, blas_column const& genos);

    double sum() const;
    double sum_col(int index) const;
    double sum_row(int index) const;
    double average() const;
};

#endif // CONTINGENCY_HPP
