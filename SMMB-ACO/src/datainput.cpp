/*
 * datainput.cpp
 *
 *  Created on: 9 nov. 2018
 *      Author: Louison Fresnais, François Courtin M2BB
 *      Project: SMMB-ACO and Genetic Algorithm for epistasis detection
 *      Under the supervision of Christine Sinoquet(Nantes University)
 *  Modified on: 19 nov 2018
 */

#include "includes/datainput.hpp"
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
	cout << "file = " << filename;
	nrows = count_rows() - header_nrows;
	ncols = count_cols();

}
/*
 * *************************************************************
 */

//Destructor
Data_input::~Data_input() {
	// TODO Auto-generated destructor stub
}
/*
 * *************************************************************
 */

/*
 *The read method iterate over the lines of the file given by the filename attribute. This is the main method of the datainput class.
 */

blas::matrix<int> Data_input::read()
{
	blas::matrix<int> matrix(nrows,ncols);
	ifstream content(filename);

	unsigned int row_pos = 1;
	unsigned int row_matrix = 0;
	//if the file is open, we iterate over its content
	if(content.is_open())
	{

		string line = "";
		while (getline(content,line))
		{
			//if the line number (first line, second line of the file for example) is in the header section, we ignore it.
			if(row_pos <= header_nrows)
			{
				cout << "HEADER Line, IGNORED \n";
			}else
			{
				unsigned int col_pos = 0;
				string temp_string = line;
				//we iterate over the line content stored in the temp_string variable
				for(auto it=temp_string.begin(); it!=temp_string.end(); ++ it)
						{
							//if the current char is not a sep, we convert it from string to int
							//and add its value in the matrix at the appropriate coordinate
							if(*it == sep){continue;}
							string temp_char(1,*it);
							int temp_value = stoi(temp_char);
							matrix(row_matrix,col_pos) = temp_value;
							col_pos++;
						}
				row_matrix++;

			}
			row_pos++;
		}

	}else
	{
		cerr << "error in file parsing \n";
		exit(-1);
	}
	content.close();
	return matrix;

}
/*
 * *************************************************************
 */

/*
 *The get_snps function will return the name of the snps in the dataset. In the final version we don't use it anymore in order because
 *of differences in snps names from naive data that are incompatible with the evaluation script.
 */

vector<string> Data_input::get_snps()
{
	vector<string> list_snps;
	ifstream content(filename);
	string temp_char;

	unsigned int row_pos = 1;
	unsigned int row_matrix = 0;
	//if the file is open, we iterate over its content.
	if(content.is_open())
	{

		string line = "";
		while (getline(content,line))
		{
			if(row_pos <= header_nrows)
			{
				unsigned int col_pos = 0;
				string temp_string = line;

				for(auto it=temp_string.begin(); it!=temp_string.end(); ++ it)
						{
							//if we arrive at a separator, it means that the previous value is a SNP, that is why we add the temp_char
							if(*it == sep)
							{
								list_snps.push_back(temp_char);
								temp_char = "";
								continue;
							}
							//add the current value of iterator to the temp_char variable.
							temp_char+= *it;
							col_pos++;
						}
				//Push back the last snp in the list_snps
				list_snps.push_back(temp_char);
				row_matrix++;
			}

			row_pos++;
		}

	}else
	{
		cerr << "error in file parsing \n";
		exit(-1);
	}
	content.close();

	return list_snps;

}
/*
 * *************************************************************
 */

/*
 *The count_rows function will count the number of rows in the file in order to remove the header when necessary
 *and initialize the storage matrix. called in the constructor
 */

unsigned int Data_input::count_rows()
{
	unsigned int nrows=0;
	string line = "";
	ifstream content(filename);
	//if the file is open, iterate over its contennt and count lines.
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
}
/*
 * *************************************************************
 */

/*
 *The count_cols function will count the number of columns in the file in order to define the storage matrix.
 *called in the constructor
 */
unsigned int Data_input::count_cols()
{
	unsigned int ncols=0;
	string line = "";
	ifstream content(filename);
	int current_line = 0;

	if(content.is_open())
	{
		while (getline(content,line))
		{
			if(current_line == 3) //We want to read only a genotype line not a header line (header is only the first line) so we pick 3rd line
			{

				string temp_string = line;

				for(auto it=temp_string.begin(); it!=temp_string.end(); ++ it)
				{
					if(*it == ','){continue;}
					ncols++;
				}
			}
			current_line++;
		}
	}else
	{
		cerr << "Error in file reading \n";
		exit(-1);
	}

	return ncols;
}
