/*
 * Miscellaneous.hpp
 *
 *  Created on: 24 nov. 2018
 *      Author: louison
 */

#ifndef MISCELLANEOUS_HPP_
#define MISCELLANEOUS_HPP_

#include <fstream>
#include <iostream>
#include <string>
#include <list>

using namespace std;

class Miscellaneous {
public:
	Miscellaneous();
	virtual ~Miscellaneous();
	void append_list(list<unsigned>, list<unsigned>);
	void remove_list_from_list(list<unsigned>, list<unsigned>);
	void append_to_file(string, string);
private:
	string filename;
	string data_to_append;
	list<unsigned> list1;
	list<unsigned> list2;

};

#endif /* MISCELLANEOUS_HPP_ */
