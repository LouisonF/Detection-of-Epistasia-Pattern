/*
 * datainput.h
 *
 *  Created on: 9 nov. 2018
 *      Author: louison
 */

#ifndef DATAINPUT_H_
#define DATAINPUT_H_

namespace std {
#include<boost>


class data_input {
public:
	//Attributes
    matrix<int> input;
	//Methods
	data_input();
	virtual ~data_input();
	file_parsing();
	
	

};

} /* namespace std */

#endif /* DATAINPUT_H_ */
