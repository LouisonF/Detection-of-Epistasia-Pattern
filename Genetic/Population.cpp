/*
 * Population.cpp
 *
 *  Created on: 19 nov. 2018
 *      Author: courtin
 */
#include <string>
#include <iostream>
#include <set>
#include "TypeDef.h"
#include "Population.h"

using namespace std;

Population::Population(int_matrix_type Mgeno, int_matrix_type Mpheno, int len_pop, int len_pattern) : InitialMatrix(Mgeno, Mpheno), nb_sol(len_pop), len_pattern(len_pattern)
{
	double_matrix_type Mpop(nb_sol, len_pattern+2); //row : sol index		col : snp 1 / snp 2 / g2 / pval
	Mpop_geno = Mpop;

}

Population::~Population() {
	// TODO Auto-generated destructor stub
}

/*
void Population::init_pop_geno(){
	int nb_snp = int(Mgeno.size2());
	vector<string> list_of_pattern (nb_sol);
	string pattern;
	int snp;


	list_of_pattern.push_back("present");
	for (int i = 0; i < nb_sol; i++){
		cout << "initialization sol n°" << i << endl;
		pattern = "present";
		while (count(list_of_pattern.begin(), list_of_pattern.end(), pattern) != 0){
			pattern.clear();
			for (int j = 0; j < len_pattern; j++){
				snp = rand()%nb_snp;  //hasard selection of n snp from an individu to form the solution
				while(count(pattern.begin(), pattern.end(), ) != 0){
					cout << "rentre" << endl;
					snp = rand()%nb_snp;
				}
				pattern.append(to_string(snp));
			}
			sort(pattern.begin(), pattern.end());
		}
		cout << pattern << endl;
		list_of_pattern.push_back(pattern);
		for (int k = 0; k < len_pattern; k++){
			Mpop_geno (i, k) = int(pattern[k]-'0'); //convert string into int
			Mpop_geno (i,len_pattern) = 0;
			Mpop_geno (i,len_pattern+1) = 0;
		}
	}
}

*/
void Population::init_pop_geno(){
	cout << "Population initialization..." << endl;

	int nb_snp = int(Mgeno.size2());
	vector<vector<int>> list_of_pattern(nb_sol);
	vector<int> pattern(len_pattern);
	int snp;

	list_of_pattern.push_back(pattern);
	for (int i = 0; i < nb_sol; i++){
		cout << "sol n°" << i+1 << " / ";
		pattern = {};
		while (count(list_of_pattern.begin(), list_of_pattern.end(), pattern) != 0){
			pattern.clear();
			for (int j = 0; j < len_pattern; j++){
				snp = rand()%nb_snp;  //hasard selection of n snp from an individu to form the solution
				while(count(pattern.begin(), pattern.end(), snp) != 0){
					snp = rand()%nb_snp;
				}
				pattern.push_back(snp);
			}
			sort(pattern.begin(), pattern.end());
		}
		list_of_pattern.push_back(pattern);
		for (int k = 0; k < len_pattern; k++){
			Mpop_geno (i, k) = pattern[k]; //convert string into int
			Mpop_geno (i,len_pattern) = 0;
			Mpop_geno (i,len_pattern+1) = 0;
		}
	}
}

void Population::display_geno_sol(){
	cout << endl << " population : " << Mpop_geno << endl;
}

double_matrix_type Population::get_Mpop_geno(){
	return(Mpop_geno);
}

void Population::set_Mpop_geno(int row, int col, double val){
	Mpop_geno(row, col) = val;
}
