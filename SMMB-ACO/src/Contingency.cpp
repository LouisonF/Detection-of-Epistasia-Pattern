/*
 * Contingency.cpp
 *
 *  Created on: 9 nov. 2018
 *      Author: Louison Fresnais, François Courtin
 *      Project: SMMB-ACO and Genetic Algorithm for epistasis detection
 *      Under the supervision of Christine Sinoquet(Nantes University)
 *      Imported from Clément Niel's code
 *  Modified on: 05 fev 2018
 */

//Module from Clement Niel's code.
#include "includes/Contingency.hpp"
#include <boost/numeric/ublas/matrix_proxy.hpp>

using namespace std;
namespace blas= boost::numeric::ublas;

//-----------------------------------------
// Constructors
//-----------------------------------------
Contingency::Contingency() : blas_dmatrix(2,3)
{
    for(unsigned i=0; i<size1(); ++i)
    {
        for(unsigned j=0; j<size2(); ++j)
            this->at_element(i,j) = 0;
    }
}
/*
 * *************************************************************
 */

Contingency::Contingency(int a, int b) : blas_dmatrix(a,b)
{
//    rnd::mt19937 rng(time(0));
//    rnd::uniform_int_distribution<> dist(30,100);
    for (unsigned i = 0; i < size1(); ++i)
    {
        for (unsigned j = 0; j < size2(); ++j)
            this->at_element(i, j) = 0; //dist(rng);
    }
}
/*
 * *************************************************************
 */

//Copy constructors
Contingency::Contingency(blas_matrix const& m) : blas_dmatrix(m.size1(), m.size2())
{
    for (unsigned i = 0; i < size1(); ++i)
    {
        for (unsigned j = 0; j < size2(); ++j)
            this->at_element(i, j) = m(i,j);
    }
}
/*
 * *************************************************************
 */

Contingency::Contingency(Contingency const& m) : blas_dmatrix(m.size1(), m.size2())
{
    for (unsigned i = 0; i < size1(); ++i)
    {
        for (unsigned j = 0; j < size2(); ++j)
            this->at_element(i, j) = m(i,j);
    }
}
/*
 * *************************************************************
 */

Contingency::Contingency(blas_column const& phenos, blas_column const& genos) : blas_dmatrix(2, 3)
{
    if(phenos.size() != genos.size())
    {
        cerr << "Genos & phenos aren't the same size. Exit" << endl;
        exit(-1);
    }

    // Initialisation to zeros
    for (unsigned i=0; i < size1(); ++i)
    {
        for (unsigned j=0; j < size2(); ++j)
            this->at_element(i, j) = 0;
    }

    // Fill contingency table
    for(unsigned i=0; i<phenos.size(); ++i)
    {
        int row_contingency = phenos(i);
        int col_contingency = genos(i);
        if((row_contingency != 0 && row_contingency != 1) || (col_contingency != 0 && col_contingency != 1 && col_contingency != 2))
        {
            cout << "Index error while building contingency table." << endl;
            continue;
        }
        this->at_element(row_contingency, col_contingency) += 1;
    }
}
/*
 * *************************************************************
 */

//-----------------------------------------
// Contingency : sum
//-----------------------------------------

//Used in the g2_test_indep run method. Will compute sum of contingency tables
double Contingency::sum() const
{
    double s=0;
    blas_dmatrix::const_iterator1 it1;
    blas_dmatrix::const_iterator2 it2;

    for(it1 = begin1(); it1 != end1(); ++it1)
    {
        for(it2 = it1.begin(); it2 != it1.end(); ++it2)
            s += *it2;
    }
    return s;
}
/*
 * *************************************************************
 */

//-----------------------------------------
// Contingency : sum_col
//-----------------------------------------

//Used in the g2_test_indep run method. Will compute sum of columns in a contingency table
double Contingency::sum_col(int index) const
{
    double s=0;
    blas_dmatrix::const_iterator1 it1;
    blas_dmatrix::const_iterator2 it2;

    for(it1 = begin1(); it1 != end1(); ++it1)
    {
        it2 = it1.begin() + index;
        s += *it2;
    }
    return s;
}
/*
 * *************************************************************
 */

//-----------------------------------------
// Contingency : sum_row
//-----------------------------------------

//Used in the g2_test_indep run method. Will compute sum of rows in a contigency table
double Contingency::sum_row(int index) const
{
    double s=0;
    blas_dmatrix::const_iterator1 it1 = begin1() + index;
    blas_dmatrix::const_iterator2 it2;

    for(it2 = it1.begin(); it2 != it1.end(); ++it2)
        s += *it2;

    return s;
}
