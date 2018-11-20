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


using namespace std;


int main () {
	srand (time(NULL));

	int_matrix_type Mgeno (50, 10);
	int_matrix_type Mpheno (50, 1);

	int nb_indiv_geno = Mgeno.size1();
	int nb_snp_pheno = Mgeno.size2();

	for (int i = 0; i < nb_indiv_geno; ++ i)
		for (int j = 0; j < nb_snp_pheno; ++ j)
			Mgeno (i, j) = rand()%3;

	int nb_indiv_pheno = Mpheno.size1();

	for (int i = 0; i < nb_indiv_pheno; ++ i)
		Mpheno (i, 0) = rand()%2;

	Population population(Mgeno, Mpheno, 20, 3);
	population.display_geno();
	population.display_pheno();
	population.init_pop_geno();
	population.display_geno_sol();
	Parent parents(Mgeno, Mpheno, 20, 3, 4);
	parents.parents_selection();
	parents.display_parents();
}
//Ressortir chaque solution avec son score et évetuellment leur p value
