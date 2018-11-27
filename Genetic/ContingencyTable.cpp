/*
 * ContengencyTable.cpp
 *
 *  Created on: 21 nov. 2018
 *      Author: courtin
 */

#include "ContingencyTable.h"
#include <string>

ContingencyTable::ContingencyTable(int_matrix_type Mgeno, int_matrix_type Mpheno, int_matrix_type Msol_geno, int len_pattern) : len_pattern(len_pattern), Mgeno(Mgeno), Mpheno(Mpheno), Msol_geno(Msol_geno){
	nb_col = pow(3, len_pattern);
	list_pattern = {};
	int_matrix_type table_cont_temp(2, nb_col);
	cont_table = table_cont_temp;

}

ContingencyTable::~ContingencyTable() {
	// TODO Auto-generated destructor stub
}

void ContingencyTable::set_pattern_list(){

	int x1 = 0;
	int x2;
	if (len_pattern == 2){
		while (x1 != 3){
			for (x2 = 0; x2 < 3; x2++){
				list_pattern.push_back(to_string(x1) + to_string(x2));
			}
			x1++;
		}
	}
	if(len_pattern == 3){
		int x3;
		while (x1 != 3){
			for (x2 = 0; x2 < 3; x2++){
				for (x3 = 0; x3 < 3; x3 ++){
					list_pattern.push_back(to_string(x1) + to_string(x2) + to_string(x3));
				}
			}
			x1++;
		}
	}
}

void ContingencyTable::set_table(){
	for (int i = 0; i < int(cont_table.size1()); i++){
		for (int j = 0; j < int(cont_table.size2()); j++){
			cont_table(i, j) = 0;
		}
	}
string pattern;
vector<string>::iterator it;
int index;

	for (int i = 0; i < int(Mpheno.size1()); i++){
		pattern = "";
		for (int j = 0; j < len_pattern; j++){
			pattern += to_string(int(Mgeno(i, Msol_geno(0,j))));
		}
		it = find(list_pattern.begin(), list_pattern.end(), pattern);
		index = distance(list_pattern.begin(), it);
		switch (Mpheno(i, 0)){
		case 0:
			cont_table(0, index)++;
			break;
		case 1:
			cont_table(1, index)++;
			break;
		}
	}
}

void ContingencyTable::display_table(){
	cout << "contingency table : " << cont_table << endl;
}

int_matrix_type ContingencyTable::get_cont_table(){
	return(cont_table);
}
