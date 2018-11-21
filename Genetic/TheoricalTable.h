/*
 * TheoricalTable.h
 *
 *  Created on: 21 nov. 2018
 *      Author: courtin
 */

#ifndef THEORICALTABLE_H_
#define THEORICALTABLE_H_

#include "TypeDef.h"
using namespace std;

class TheoricalTable {
public:
	int nb_col;
	int len_pattern;
	int_matrix_type cont_table;
	int_matrix_type theo_table;
public:
	TheoricalTable(int, int_matrix_type);
	virtual ~TheoricalTable();
	float sum_row(int);
	float sum_col(int);
	void set_table();
	void display_table();
};

#endif /* THEORICALTABLE_H_ */
