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
	Parent parents(Mgeno, Mpheno, len_pop, len_pattern, nb_parents);
	parents.parents_selection();
	parents.display_parents();

	int_matrix_type Msol_geno(1, len_pattern);


	for (int i = 0; i < len_pop; i++){
		for (int j = 0; j < len_pattern; j++){
			Msol_geno(0, j) = population.Mpop_geno(i, j);
		}
		ContingencyTable cont_table(Mgeno, Mpheno, Msol_geno, len_pattern);
		cont_table.set_pattern_list();
		cont_table.set_table();
		cont_table.display_table();

		TheoricalTable theo_table(len_pattern, cont_table.cont_table);
		theo_table.set_table();
		theo_table.display_table();
	}
}
//Ressortir chaque solution avec son score et évetuellment leur p value
