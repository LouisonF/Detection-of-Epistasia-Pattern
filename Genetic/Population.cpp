/*
 * Population.cpp
 *
 *  Created on: 19 nov. 2018
 *      Author: courtin
 */

#include "Population.h"
using namespace std;

Population::Population(matrix_of_int_matrix_type Mgeno, int_matrix_type Mpheno, int len_pop, int len_patern) {
	nb_sol = len_pop;
	nb_indiv = Mgeno.size1();
	select_sol(nb_sol, 1);

	for (int i = 0; i < nb_sol; i++){
		select_sol(i, 0) = rand()%((nb_indiv));
	}

	int nb_snp = Mgeno.size2();

	int_matrix_type Msol (nb_indiv, len_patern); //Matrix of a solution

	Mpop_geno(nb_sol, 1);

	int nb_indiv_pop = Msol.size1();

    for (int i = 0; i < nb_sol; i++){
    	for (int j = 0; j < nb_indiv_pop; j++){
    		for (int k = 0; k < len_patern; k++){
            	Msol (j, k) = Mgeno(select_sol(i, 0), rand()%(nb_snp)); //hasard selection of 3 snp from an individu to forme the solution
            	Msol (j, k) = Mgeno(select_sol(i, 0), rand()%(nb_snp));
            	Msol (j, k) = Mgeno(select_sol(i, 0), rand()%(nb_snp));
    		}
    	}
    	Mpop_geno(i, 0) = Msol; //add the solution matrix in the list of solutions
    }

    Mpop_pheno(nb_indiv_pop, 1);

    for (int i = 0; i < nb_indiv_pop; i++){
    	Mpop_pheno(i, 0) = Mpheno(select_sol(i, 0), 0);
    }

}

Population::~Population() {
	// TODO Auto-generated destructor stub
}


void Population::display_geno_sol(){
	for (int i = 0; i < nb_indiv; i++)
		cout << " genno solution " << i+1 << ":" << Mpop_geno(i, 0) << endl;
}
void Population::display_pheno_sol(){
	for (int i = 0; i < nb_indiv; i++)
		cout << " pheno solution " << i+1 << ":" << Mpop_pheno(i, 0) << endl;
}
