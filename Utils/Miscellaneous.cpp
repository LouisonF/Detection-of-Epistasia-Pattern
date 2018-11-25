/*
 * Miscellaneous.cpp
 *
 *  Created on: 24 nov. 2018
 *      Author: louison
 */

#include "Miscellaneous.hpp"
using namespace std;

Miscellaneous::Miscellaneous() {
	// TODO Auto-generated constructor stub

}

Miscellaneous::~Miscellaneous() {
	// TODO Auto-generated destructor stub
}

void Miscellaneous::append_list(list<unsigned> list1, list<unsigned> list2)
{
	list<unsigned> buffer = list2;
	list1.insert(list1.end(),buffer.begin(), buffer.end());
}

void Miscellaneous::remove_list_from_list(list<unsigned> list1, list<unsigned> list2)
{
	for(auto it=list1.begin(); it!=list1.end(); ++ it)
	{
		unsigned temp = *it;
		list2.remove(temp);
	}
}

void Miscellaneous::append_to_file(string filename, string data_to_append)
{
	ofstream content;
	content.open(filename, ofstream::out | ofstream::app);
	//option out allow to open the file in writing mode.
	//option app allow to append data at the end of the file.

	if(content.is_open())
	{
		content << data_to_append;
	}
	else
	{
		cerr << "error in file parsing \n";
		exit(-1);
	}
	content.close();
}
