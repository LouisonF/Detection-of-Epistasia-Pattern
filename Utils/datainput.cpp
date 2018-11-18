/*
 * datainput.cpp
 *
 *  Created on: 9 nov. 2018
 *      Author: louison
 */

#include "datainput.hpp"
#include <stdlib.h>

#include <iostream>
#include <string>
#include <fstream>
#include <algorithm>
#include <map>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>



using namespace std;

//Constructor

Data_input::Data_input(const string filename, char sep, const unsigned int header_nrows) : filename(filename), sep(sep), header_nrows(header_nrows)
{
	nrows = count_rows() - header_nrows;
	ncols = count_cols();

	//bool transpose = false; Temporary, removed when i'll decide if i add the transpose mode or not.


}
//Destructor
Data_input::~Data_input() {
	// TODO Auto-generated destructor stub
}
//Debut des méthodes
//***********************************************************************************************************

 //This function read the input file and store data in a matrix.

void Data_input::read()
{
	blas::matrix<int> matrix(nrows,ncols);
	ifstream content(filename);

	unsigned int row_pos = 1;
	unsigned int row_matrix = 0;

	if(content.is_open())
	{

		string line = "";
		while (getline(content,line))

			if(row_pos <= header_nrows)
			{
				cout << "HEADER Line, IGNORED \n";
			}else
			{
				unsigned int col_pos = 0;
				string temp_string = line;

				for(auto it=temp_string.begin(); it!=temp_string.end(); ++ it)
						{
							if(*it == sep){continue;}
							string temp_char(1,*it);
							int temp_value = stoi(temp_char);
							matrix(row_matrix,col_pos) = temp_value;
							col_pos++;
						}
				row_matrix++;

			}
			row_pos++;


	}else
	{
		cerr << "error in file parsing \n";
		exit(-1);
	}
	cout << matrix <<endl;
	content.close();
}

// This function count the number of rows in the input file
unsigned int Data_input::count_rows()
{
	unsigned int nrows=0;
	string line = "";
	ifstream content(filename);

	if(content.is_open())
	{
		while(getline(content,line))
		{
			nrows++;
		}
	}else
	{
		cerr << "Error in file reading \n";
		exit(-1);
	}
	return nrows;
};

//This function count the number of cols in the input file
unsigned int Data_input::count_cols()
{
	unsigned int ncols=0;
	string line = "";
	ifstream content(filename);

	if(content.is_open())
	{
		getline(content,line);

		string temp_string = line;

		for(auto it=temp_string.begin(); it!=temp_string.end(); ++ it)
				{
					if(*it == ','){continue;}
					ncols++;
				}

	}else
	{
		cerr << "Error in file reading \n";
		exit(-1);
	}
	return ncols;
};
