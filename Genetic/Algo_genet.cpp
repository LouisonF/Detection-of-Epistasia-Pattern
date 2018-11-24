//courtin francois
//fresnais louison
//Projet algo genetic

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/array.hpp>
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <vector>
#include "InitialMatrix.h"
#include "Population.h"
#include "Parent.h"
#include "ContingencyTable.h"
#include "TheoricalTable.h"
#include "G2test.h"


using namespace std;


int main () {
	srand (time(NULL));

	int_matrix_type Mgeno (1000, 10);
	int_matrix_type Mpheno (1000, 1);

	int nb_indiv_geno = Mgeno.size1();
	int nb_snp_pheno = Mgeno.size2();

	for (int i = 0; i < nb_indiv_geno; ++ i)
		for (int j = 0; j < nb_snp_pheno; ++ j)
			Mgeno (i, j) = rand()%3;

	int nb_indiv_pheno = Mpheno.size1();

	for (int i = 0; i < nb_indiv_pheno; ++ i)
		Mpheno (i, 0) = rand()%2;


	int len_pop = 20;
	int len_pattern = 2;
	int nb_parents = 4;

	Population population(Mgeno, Mpheno, len_pop, len_pattern);
	population.display_geno();
	population.display_pheno();
	population.init_pop_geno();
	population.display_geno_sol();

	int_matrix_type Msol_geno(1, len_pattern); //temp solution for the loop
	vector<float> G2_res;
	vector<float> pval_res;

	for (int i = 0; i < len_pop; i++){
		for (int j = 0; j < len_pattern; j++){
			Msol_geno(0, j) = population.get_Mpop_geno()(i, j);
		}
		ContingencyTable cont_table(Mgeno, Mpheno, Msol_geno, len_pattern);
		cont_table.set_pattern_list();
		cont_table.set_table();
		cont_table.display_table();

		TheoricalTable theo_table(len_pattern, cont_table.get_cont_table());
		theo_table.set_table();
		theo_table.display_table();

		G2test G2(cont_table.get_cont_table(), theo_table.get_theo_table());
		G2.run_G2();
		G2.display_g2();
		G2_res.push_back(G2.get_g2());
		pval_res.push_back(G2.get_pval());
	}
	vector<float> oui = sort(G2_res.begin(), G2_res.end());
	Parent parents(Mgeno, Mpheno, len_pop, len_pattern, nb_parents);
	parents.parents_selection();
	parents.display_parents();
}
//Calculer ma médiane depuis le vecteur trié, faire la sélection des parents sur des valeurs uniquement au dessus de cette médiane.
//Verifier si sort trie bien en croissant.
