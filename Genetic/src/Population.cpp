/*
 * Population.cpp
 *
 *  Created on: 19 nov. 2018
 *      Author: Louison Fresnais, François Courtin
 *      Project: SMMB-ACO and Genetic Algorithm for epistasis detection
 *      Under the supervision of Christine Sinoquet(Nantes University)
 */
#include <string>
#include <iostream>
#include <set>
#include "TypeDef.h"
#include "Population.h"

using namespace std;

Population::Population(int_matrix_type Mgeno, int_matrix_type Mpheno, int &len_pop, int len_pattern) : Mgeno(Mgeno), Mpheno(Mpheno), nb_sol(len_pop), len_pattern(len_pattern)
{
	double_matrix_type Mpop(nb_sol, len_pattern+2); //row : sol index		col : snp 1 / snp 2 / g2 / pval
	Mpop_geno = Mpop;

}

Population::~Population() {
	// TODO Auto-generated destructor stub
}

/*
 * *************************************************************
 */
//recursive factorial function of number "N"
long double Population::fact(long int N){
	if(N<=1)
	return 1;
	else
	return(N*fact(N-1));
}

/*
 * *************************************************************
 */

//Check the max number of unique solutions that can be put in the initial population
void Population::check_pop_len(){
	long double max_len;
	max_len = fact(Mgeno.size2())/(fact(Mgeno.size2()-len_pattern)*fact(len_pattern));
	//If the length of the population is higher than the max solutions then set the length to the max solutions
	if (max_len < nb_sol)
	{
		cout << "Le nombre de solutions demandées dans la population est trop grand." << endl << "Le nombre maximal de solutions possible est " << max_len << endl;
		cout <<"La taille de la population est fixée à " << max_len << endl;
		nb_sol = max_len;
	}
}

/*
 * *************************************************************
 */

//Initialize the population with uniques solutions
void Population::init_pop_geno(){
	//cout << "Population initialization..." << endl;

	int nb_snp = int(Mgeno.size2());
	vector<vector<int>> list_of_pattern(nb_sol);
	vector<int> pattern(len_pattern);
	int snp;

	list_of_pattern.push_back(pattern);
	for (int i = 0; i < nb_sol; i++){
		pattern = {};
		//Pick a solution while it isn't already in the populations
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
			Mpop_geno (i, k) = pattern[k];
			Mpop_geno (i,len_pattern) = 0; //Set the g2 score associated to the solution to 0
			Mpop_geno (i,len_pattern+1) = 1; //Set the p-value associated to the solution to 1
		}
	}
}

/*
 * *************************************************************
 */

void Population::display_geno_sol(){
	cout << endl << " population : " << Mpop_geno << endl;
}

/*
 * *************************************************************
 */

double_matrix_type Population::get_Mpop_geno(){
	return(Mpop_geno);
}

/*
 * *************************************************************
 */

//Set new values to the population
void Population::set_Mpop_geno(int row, int col, double val){
	Mpop_geno(row, col) = val;
}
