/*
 * TheoricalTable.cpp
 *
 *  Created on: 21 nov. 2018
 *      Author: Louison Fresnais, François Courtin
 *      Project: SMMB-ACO and Genetic Algorithm for epistasis detection
 *      Under the supervision of Christine Sinoquet(Nantes University)
 */

#include "includes/TheoricalTable.h"

TheoricalTable::TheoricalTable(int len_pattern, int_matrix_type cont_table) : len_pattern(len_pattern), cont_table(cont_table){
	nb_col = pow(3, len_pattern);
	int_matrix_type theo_table_temp(2, nb_col);
	theo_table = theo_table_temp;

}

TheoricalTable::~TheoricalTable() {
	// TODO Auto-generated destructor stub
}

/*
 * *************************************************************
 */

//Function to compute the sum of a row
double TheoricalTable::sum_row(int row){
	double sum = 0;
	for (int i = 0; i < nb_col; i++){
		sum += cont_table(row, i);
	}
	return(sum);
}

/*
 * *************************************************************
 */

//Function to compute the sum of a col
double TheoricalTable::sum_col(int col){
	double sum = 0;
	for (int i = 0; i < int(cont_table.size1()); i++){
		sum += cont_table(i, col);
	}
	return(sum);
}

/*
 * *************************************************************
 */

//Set the theorical contingency table based on the contingency table
void TheoricalTable::set_table(){
	//cout << "Set theorical table..." << endl;
	double total = sum_row(0) + sum_row(1);
	for (int i = 0; i < int(theo_table.size1()); i++){
		for (int j = 0; j < int(theo_table.size2()); j++){
			theo_table(i, j) = (sum_row(i) * sum_col(j))/total;
		}
	}
}

/*
 * *************************************************************
 */

void TheoricalTable::display_table(){
	cout << "Theorical table : " << theo_table << endl;
}

/*
 * *************************************************************
 */

int_matrix_type TheoricalTable::get_theo_table(){
	return(theo_table);
}
