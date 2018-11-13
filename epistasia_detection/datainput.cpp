/*
 * datainput.cpp
 *
 *  Created on: 9 nov. 2018
 *      Author: louison
 */

#include "datainput.hpp"

#include <iostream>
#include <string>
#include <fstream>
#include <map>
#include <boost/numeric/ublas/matrix.hpp>



using namespace std;

//Constructor
/*
Data_input::Data_input(const string & filename, char sep, const unsigned int header_nrows, bool transpose)
{
	nrows = nrows = count_rows() - header_nrows;
	ncols = count_cols();

	if(!transpose)
		read();
	else
		read_transpose();
}*/
//Temporary constructor
Data_input::Data_input() {
	// TODO Auto-generated constructor stub

}

Data_input::~Data_input() {
	// TODO Auto-generated destructor stub
}
//Debut des méthodes
//***********************************************************************************************************
/*Data_input::data() //Not sure about the utility of this, soon deprecated
{
	return this.matrix;
}*/
 //File reading

void Data_input::read(unsigned int header_nrows,unsigned int nrows,unsigned intncols)
{
	blas::matrix<int> matrix(nrows,ncols);
	ifstream file_content("test_data.txt");

	unsigned int row_pos = 1;

	if(file_content.is_open())
	{

		std::string line = "";
		while (std::getline(file_content,line)) //Work but don't display in eclipse consol ... check with terminal
		{
			if(row_pos <= header_nrows)
			{
				cout << "HEADER Line, IGNORED \n";
			}else
			{
				cout << line << endl;
				unsigned int col_pos = 1;
				for(std::string::iterator it = line.begin(); it != line.end();it++) // Cette boucle ne fonctionne pas
				{
					//matrix (row_pos,col_pos) = *it;
					col_pos++;
					//test
					cout << "char" + *it << endl;
					//cout << matrix(row_pos,col_pos);
				}
			}
			row_pos++;
			cout << row_pos << endl;

		}
file_content.close();
	}else
	{
		cerr << "error in file parsing \n";
		exit(1);
	}
};
/*
unsigned int Data_input::count_rows()
{
	unsigned int nrows=0;
	string line;
	ifstream file_content(filename);

	if(file)
	{
		while(getline(file_content,line))
		{
			nrows++;
		}
	}else
	{
		cerr << "Error in file reading \n";
	}
	return nrows;
}
*/
