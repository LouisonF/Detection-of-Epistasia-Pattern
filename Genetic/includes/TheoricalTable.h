/*
 * TheoricalTable.h
 *
 *  Created on: 21 nov. 2018
 *      Author: Louison Fresnais, François Courtin
 *      Project: SMMB-ACO and Genetic Algorithm for epistasis detection
 *      Under the supervision of Christine Sinoquet(Nantes University)
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
	double sum_row(int);
	double sum_col(int);
	void set_table();
	void display_table();
	int_matrix_type get_theo_table();
};

#endif /* THEORICALTABLE_H_ */
