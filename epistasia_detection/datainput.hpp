/*
 * datainput.h
 *
 *  Created on: 9 nov. 2018
 *      Author: louison
 */

#ifndef DATAINPUT_HPP_
#define DATAINPUT_HPP_

using namespace std;
#include <iostream>
#include <string>
#include <fstream>
#include <map>
#include <boost/numeric/ublas/matrix.hpp>


namespace blas=boost::numeric::ublas;

class Data_input {
public:
Data_input(const string filename, char sep, const unsigned int header_nrows);
~Data_input();
    /*!
     * \brief Get the imported matrix
     * \return Imported Matrix
     */
    blas::matrix<int> data();

    /*!
     * \brief Count rows in the CSV file
     * \return Row count
     */
    unsigned int count_rows( );

    /*!
     * \brief Count cols in the CSV file
     * \return Column count
     */
    unsigned int count_cols( );

    /*!
     * \brief Import the CSV data in _matrix (Matrix<T> class)
     */
    void read(); //A changer

    /*!
     * \brief Transpose and import the CSV data in _matrix (Matrix<T> class)
     */
    void read_transpose();


private:
    string filename;
    char sep;
    unsigned int header_nrows;
    unsigned int nrows;
    unsigned int ncols;
    blas::matrix<int> matrix();
};


#endif /* DATAINPUT_HPP_ */
