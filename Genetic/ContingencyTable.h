/*
 * ContengencyTable.h
 *
 *  Created on: 21 nov. 2018
 *      Author: courtin
 */

#ifndef CONTINGENCYTABLE_H_
#define CONTINGENCYTABLE_H_
#include "TypeDef.h"

using namespace std;

class ContingencyTable {
private:
	int len_pattern;
	int nb_col;
	vector<vector<int>> list_pattern;
	vector<int> pattern;
	int_matrix_type Mgeno;
	int_matrix_type Mpheno;
	int_matrix_type Msol_geno;
	int_matrix_type cont_table;
public:
	ContingencyTable(int_matrix_type, int_matrix_type, int_matrix_type, int);
	virtual ~ContingencyTable();
	void set_pattern_list();
	void set_table();
	void display_table();
	int_matrix_type get_cont_table();
};

#endif /* CONTINGENCYTABLE_H_ */
