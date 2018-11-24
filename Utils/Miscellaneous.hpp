/*
 * Miscellaneous.hpp
 *
 *  Created on: 24 nov. 2018
 *      Author: louison
 */

#ifndef MISCELLANEOUS_HPP_
#define MISCELLANEOUS_HPP_

#include <fstream>

class Miscellaneous {
public:
	Miscellaneous();
	virtual ~Miscellaneous();
	void append_list(list<unsigned>, list<unsigned>);
	void remove_list_from_list(list<unsigned>, list<unsigned>);
	void append_to_file(string, string)
};

#endif /* MISCELLANEOUS_HPP_ */
