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
	int_matrix_type Mpop(nb_sol, len_pattern);
	Mpop_geno = Mpop;

}

Population::~Population() {
	// TODO Auto-generated destructor stub
}


void Population::init_pop_geno(){
	int nb_snp = int(Mgeno.size2());
	set<string> list_of_pattern = {"present"};
	string pattern;
	int snp;

    for (int i = 0; i < nb_sol; i++){
    	pattern = "present";
    	while (list_of_pattern.count(pattern) != 0){
    		pattern.clear();
    		for (int j = 0; j < len_pattern; j++){
    			snp = rand()%nb_snp;  //hasard selection of n snp from an individu to form the solution
    			pattern.append(to_string(snp));
    		}
    	}
    	list_of_pattern.emplace(pattern);
		for (int k = 0; k < len_pattern; k++){
			Mpop_geno (i, k) = int(pattern[k]-'0'); //convert string into int
		}
    }
}

void Population::display_geno_sol(){
	cout << " population : " << Mpop_geno << endl;
}

int_matrix_type Population::get_Mpop_geno(){
	return(Mpop_geno);
}

