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
	ifstream content("test_data.txt");
	char c;

	unsigned int row_pos = 1;

	if(content.is_open())
	{

		std::string line = "";
		while (std::getline(content,line)) //Work but don't display in eclipse consol ... check with terminal
		{
			//line = line.erase(',');
			//if(row_pos <= header_nrows)
			//{
			//	cout << "HEADER Line, IGNORED \n";
			//}else

				//cout << line << endl;
				unsigned int col_pos = 1;
				string temp_string = line;
				//for(std::string::size_type i = 0; i < ss.length();i++) // Cette boucle ne fonctionne pas
				for(auto it=temp_string.begin(); it!=temp_string.end(); ++ it, ++col_pos)
						{
							if(*it == ','){continue;}
							cout << "char is ";
							//cout <<  *it << endl;
							string temp_char(1,*it);
							cout << stoi(temp_char) << endl;
							//matrix (row_pos,col_pos) = std::stoi(temp_char); // TODO: enter value in matrix, not working atm
						}
			row_pos++;
			cout << row_pos << endl;

		}
content.close();
	}else
	{
		cerr << "error in file parsing \n";
		exit(-1);
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
