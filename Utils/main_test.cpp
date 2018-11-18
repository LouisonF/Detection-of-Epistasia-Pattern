/*
 * main_test.cpp
 *
 *  Created on: 13 nov. 2018
 *      Author: louison
 */
#include "datainput.hpp"
#include "Parametersfileparsing.hpp"

#include <iostream>
#include <string>
#include <fstream>
#include <map>
#include <boost/numeric/ublas/matrix.hpp>

using namespace std;

int main()
{
	string filename = "test_data.txt";
	unsigned int header_nrows = 2;
	char sep = ',';
	Data_input test(filename, sep, header_nrows);
	test.read();

}

