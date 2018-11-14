/*
 * main_test.cpp
 *
 *  Created on: 13 nov. 2018
 *      Author: louison
 */
#include "datainput.hpp"

#include <iostream>
#include <string>
#include <fstream>
#include <map>
#include <boost/numeric/ublas/matrix.hpp>

using namespace std;

int main()
{
	unsigned int header_nrows = 2;
	unsigned int nrows = 100;
	unsigned int ncols = 100;

	Data_input test;
	test.read(header_nrows,nrows,ncols);
}

