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
void read();
unsigned int count_rows();
unsigned int count_cols();

private:
    string filename;
    char sep;
    unsigned int header_nrows;
    unsigned int nrows;
    unsigned int ncols;
    blas::matrix<int> matrix();
};


#endif /* DATAINPUT_HPP_ */
