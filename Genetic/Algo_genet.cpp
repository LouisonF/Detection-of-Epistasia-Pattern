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
#include "Child.h"


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
	int nb_parents = 2;

	Population population(Mgeno, Mpheno, len_pop, len_pattern);
	population.display_geno();
	population.display_pheno();
	population.init_pop_geno();
	population.display_geno_sol();

	int_matrix_type Msol_geno(1, len_pattern); //temp solution for the loop
	vector<float> G2_res;

	for (int i = 0; i < len_pop; i++){
		for (int j = 0; j < len_pattern; j++){
			Msol_geno(0, j) = population.get_Mpop_geno()(i, j);
		}
		ContingencyTable cont_table_pop(Mgeno, Mpheno, Msol_geno, len_pattern);
		cont_table_pop.set_pattern_list();
		cont_table_pop.set_table();
		//cont_table.display_table();

		TheoricalTable theo_table_pop(len_pattern, cont_table_pop.get_cont_table());
		theo_table_pop.set_table();
		//theo_table.display_table();

		G2test G2_pop(cont_table_pop.get_cont_table(), theo_table_pop.get_theo_table());
		G2_pop.run_G2();
		G2_pop.display_g2();
		G2_res.push_back(G2_pop.get_g2());
		population.set_Mpop_geno(i,len_pattern,G2_pop.get_g2());
		population.set_Mpop_geno(i,len_pattern+1,G2_pop.get_pval());
	}
	population.display_geno_sol();


	float median;

	sort(G2_res.begin(), G2_res.end());

	median = G2_res[G2_res.size()/2]; //Median of the solutions' G2

	////////////////////////
	//Entrée dans le while//
	////////////////////////


	Parent parents(len_pop, nb_parents, len_pattern, median, population.get_Mpop_geno());
	parents.parents_selection();
	parents.display_parents();

	Child children(population.get_Mpop_geno(), parents.get_MParents(), len_pattern);
	children.set_children();
	children.mutation();
	children.display_children();

	int_matrix_type Mchild(1, len_pattern); //Temp matrix for the loop

	for (int i = 0; i < nb_parents; i++){
		for (int j = 0; j < len_pattern; j++){
			Mchild(0, j) = children.get_MChildren()(i, j);
		}

		ContingencyTable cont_table_child(Mgeno, Mpheno, Mchild, len_pattern);
		cont_table_child.set_pattern_list();
		cont_table_child.set_table();
		cont_table_child.display_table();

		TheoricalTable theo_table_child(len_pattern, cont_table_child.get_cont_table());
		theo_table_child.set_table();
		theo_table_child.display_table();

		G2test G2_child(cont_table_child.get_cont_table(), theo_table_child.get_theo_table());
		G2_child.run_G2();
		G2_child.display_g2();
	}

}
//Comparer G2 children et parents, vérifier pval, remplacer les parents dans Mpop_geno refaire un vector des g2 recalculer la médiane et boucler.
//COMMENTER LE CODE GROS GUIGNOLE
