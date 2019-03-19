/*
 * ContengencyTable.cpp
 *
 *  Created on: 21 nov. 2018
 *      Author: Louison Fresnais, François Courtin
 *      Project: SMMB-ACO and Genetic Algorithm for epistasis detection
 *      Under the supervision of Christine Sinoquet(Nantes University)
 */

#include "includes/ContingencyTable.h"
#include <string>

ContingencyTable::ContingencyTable(int_matrix_type &Mgeno, int_matrix_type &Mpheno, int_matrix_type Msol_geno, int len_pattern) : len_pattern(len_pattern), Mgeno(Mgeno), Mpheno(Mpheno), Msol_geno(Msol_geno){
	nb_col = pow(3, len_pattern);
	int_matrix_type table_cont_temp(2, nb_col);
	cont_table = table_cont_temp;

}

ContingencyTable::~ContingencyTable() {
	// TODO Auto-generated destructor stub
}

/*
 * *************************************************************
 */

//Set a list of all possibble patterns of snp according to the length of the pattern
void ContingencyTable::set_pattern_list(){
	//cout << "Set pattern list..." << endl;

	int x1 = 0;
	int x2;
	if (len_pattern == 2){
		while (x1 != 3){
			for (x2 = 0; x2 < 3; x2++){
				pattern.clear();
				pattern.push_back(x1);
				pattern.push_back(x2);
				list_pattern.push_back(pattern);
			}
			x1++;
		}
	}
	if(len_pattern == 3){
		int x3;
		while (x1 != 3){
			for (x2 = 0; x2 < 3; x2++){
				for (x3 = 0; x3 < 3; x3 ++){
					pattern.clear();
					pattern.push_back(x1);
					pattern.push_back(x2);
					pattern.push_back(x3);
					list_pattern.push_back(pattern);
				}
			}
			x1++;
		}
	}
}

/*
 * *************************************************************
 */

//Set the contingency table
void ContingencyTable::set_table(){
	//cout << "Set contingency table..." << endl;

	//Initalize the table with 0 value everywhere
	for (int i = 0; i < int(cont_table.size1()); i++){
		for (int j = 0; j < int(cont_table.size2()); j++){
			cont_table(i, j) = 0;
		}
	}


	vector<vector<int>>::iterator it;
	int index;

	for (int i = 0; i < int(Mpheno.size1()); i++){
		pattern.clear();
		//Get the pattern of the solution i and put it into a vector
		for (int j = 0; j < len_pattern; j++){
			pattern.push_back(int(Mgeno(i, Msol_geno(0,j))));
		}
		//Get the index of the selected pattern into the list of all possible patterns
		it = find(list_pattern.begin(), list_pattern.end(), pattern);
		index = distance(list_pattern.begin(), it);
		//Incrementation of the case corresponding to the genotype of the pattern at the right index depending to the penotype
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

/*
 * *************************************************************
 */

void ContingencyTable::display_table(){
	//cout << "contingency table : " << cont_table << endl;
}

/*
 * *************************************************************
 */

int_matrix_type ContingencyTable::get_cont_table(){
	return(cont_table);
}
