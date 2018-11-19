/*
 * Population.cpp
 *
 *  Created on: 19 nov. 2018
 *      Author: courtin
 */
#include "TypeDef.h"
#include "Population.h"
using namespace std;

Population::Population(int_matrix_type Mgeno, int_matrix_type Mpheno, int len_pop, int len_patern) : Mgeno(Mgeno), Mpheno(Mpheno), nb_sol(len_pop), len_patern(len_patern)
{

	nb_indiv = Mgeno.size1();
	sol_selection();
	init_pop_geno();
	init_pop_pheno();

}

Population::~Population() {
	// TODO Auto-generated destructor stub
}

void Population::sol_selection(){
	select_sol(nb_sol, 1);
	for (int i = 0; i < nb_sol; i++){
		select_sol(i, 0) = rand()%((nb_indiv));
	}
}

void Population::init_pop_geno(){
	int nb_snp = Mgeno.size2();

	int_matrix_type Msol_geno (nb_indiv, len_patern); //Matrix of a solution

	Mpop_geno(nb_sol, 1);

	int nb_indiv_pop = Msol_geno.size1();

	int snp1;
	int snp2;
	int snp3;

    for (int i = 0; i < nb_sol; i++){
    	for (int j = 0; j < nb_indiv_pop; j++){
    		for (int k = 0; k < len_patern; k++){
    			Msol_geno (j, k) = Mgeno(select_sol(i, 0), (snp1 = rand()%(nb_snp))); //hasard selection of 3 snp from an individu to form the solution
    			Msol_geno (j, k) = Mgeno(select_sol(i, 0), (snp2 = rand()%(nb_snp)));
    			Msol_geno (j, k) = Mgeno(select_sol(i, 0), (snp3 = rand()%(nb_snp)));
    		}
    	}
    	Mpop_geno(i, 0) = Msol_geno; //add the solution matrix in the list of solutions
    	//Mpop_snp(i, 0);
    }
}

void Population::init_pop_pheno(){


	int_matrix_type Msol_pheno (nb_indiv, 1); //Matrix of a solution

	int nb_indiv_pop = Msol_pheno.size1();

    Mpop_pheno(nb_indiv_pop, 1);

	matrix_of_int_matrix_type sol_list(nb_sol, 1);

    for (int i = 0; i < nb_sol; i++){
    	for (int j = 0; j < nb_indiv_pop; j++){
    		Msol_pheno (j, 0) = Mpheno(j, 0);
    	}
    	Mpop_pheno(i, 0) = Msol_pheno; //add the solution matrix in the list of solutions
    }
}

void Population::display_geno_sol(){
	for (int i = 0; i < nb_indiv; i++)
		cout << " genno solution " << i+1 << ":" << Mpop_geno(i, 0) << endl;
}
void Population::display_pheno_sol(){
	for (int i = 0; i < nb_indiv; i++)
		cout << " pheno solution " << i+1 << ":" << Mpop_pheno(i, 0) << endl;
}
