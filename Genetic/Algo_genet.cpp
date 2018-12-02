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
#include <bits/stdc++.h>
#include <datainput.hpp>


using namespace std;


int main () {
	srand (time(NULL));

	int_matrix_type Mgeno (5000, 10);
	int_matrix_type Mpheno (5000, 1);

	int nb_indiv_geno = Mgeno.size1();
	int nb_snp = Mgeno.size2();

	for (int i = 0; i < nb_indiv_geno; ++ i)
		for (int j = 0; j < nb_snp; ++ j)
			Mgeno (i, j) = rand()%3;

	int nb_indiv_pheno = Mpheno.size1();

	for (int i = 0; i < nb_indiv_pheno; ++ i)
		Mpheno (i, 0) = rand()%2;


	int len_pop = 120;
	int len_pattern = 3;
	int nb_parents = 4;
	float alpha = 0.05;
	int nb_it = 100;
	float P_mutation = 10;
	float P_selection = 10;

	Population population(Mgeno, Mpheno, len_pop, len_pattern);
	//population.display_geno();
	//population.display_pheno();
	population.init_pop_geno();
	//population.display_geno_sol();

	int_matrix_type Msol_geno(1, len_pattern); //temp solution for the loop
	vector<float> G2_res;

	for (int i = 0; i < len_pop; i++){
		for (int j = 0; j < len_pattern; j++){
			Msol_geno(0, j) = population.get_Mpop_geno()(i, j);
		}
		ContingencyTable cont_table_pop(Mgeno, Mpheno, Msol_geno, len_pattern);
		cont_table_pop.set_pattern_list();
		cont_table_pop.set_table();
		//cont_table_pop.display_table();

		TheoricalTable theo_table_pop(len_pattern, cont_table_pop.get_cont_table());
		theo_table_pop.set_table();
		//theo_table_pop.display_table();

		G2test G2_pop(cont_table_pop.get_cont_table(), theo_table_pop.get_theo_table());
		G2_pop.run_G2();
		//G2_pop.display_g2();
		G2_res.push_back(G2_pop.get_g2());
		population.set_Mpop_geno(i,len_pattern,G2_pop.get_g2());
		population.set_Mpop_geno(i,len_pattern+1,G2_pop.get_pval());

	}
	//population.display_geno_sol();


	float median;

	sort(G2_res.begin(), G2_res.end());

	median = G2_res[G2_res.size()/2]; //Median of the solutions' G2

	////////////////////////
	//Entrée dans le while//
	////////////////////////
	int iterator = 0;
	while (iterator < nb_it){

		Parent parents(len_pop, nb_parents, len_pattern, median, P_selection, population.get_Mpop_geno());
		parents.parents_selection();
		//parents.display_parents();

		Child children(population.get_Mpop_geno(), parents.get_MParents(), len_pattern, P_mutation, nb_snp);
		children.set_children();
		children.mutation();
		//children.display_children();

		int_matrix_type Mchild(1, len_pattern); //Temp matrix for the loop

		for (int i = 0; i < nb_parents; i++){
			for (int j = 0; j < len_pattern; j++){
				Mchild(0, j) = children.get_MChildren()(i, j);
			}

			ContingencyTable cont_table_child(Mgeno, Mpheno, Mchild, len_pattern);
			cont_table_child.set_pattern_list();
			cont_table_child.set_table();
			//cont_table_child.display_table();

			TheoricalTable theo_table_child(len_pattern, cont_table_child.get_cont_table());
			theo_table_child.set_table();
			//theo_table_child.display_table();

			G2test G2_child(cont_table_child.get_cont_table(), theo_table_child.get_theo_table());
			G2_child.run_G2();
			//G2_child.display_g2();

			//COMPARAISON AVEC PARENT ET REMPLACEMENT DANS LA POPULATION
			if ((G2_child.get_g2() > population.get_Mpop_geno()(parents.get_MParents()(i,0), len_pattern)) and (G2_child.get_pval() < alpha)){
				for (int j = 0; j < len_pattern; j++){
					population.set_Mpop_geno( parents.get_MParents()(i,0) , j , children.get_MChildren()(i,j) );
				}
				population.set_Mpop_geno( parents.get_MParents()(i,0) , len_pattern , G2_child.get_g2() );
				population.set_Mpop_geno( parents.get_MParents()(i,0) , len_pattern+1 , G2_child.get_pval() );
			}
		}
		G2_res.clear();
		for (int i = 0; i < len_pop; i++){
			G2_res.push_back(population.get_Mpop_geno()(i,len_pattern));
		}
		sort(G2_res.begin(), G2_res.end());
		median = G2_res[G2_res.size()/2]; //Median of the solutions' G2

		iterator++;
		cout << "iteration n°" << iterator << endl;
	}

	//////////////
	//FIN WHILE//
	////////////

	vector<string> list_pattern;
	string pattern;
	for (int i = 0; i < len_pop; i++){
		pattern = "";
		for (int j = 0; j < len_pattern; j++){
			pattern += to_string(int(population.get_Mpop_geno()(i,j)));
		}
		if (count(list_pattern.begin(), list_pattern.end(), pattern) == 0){
			list_pattern.push_back(pattern);
		}
	}

	vector<string> list_sol;
	for (int i = 0; i < len_pop; i++){
		pattern = "";
		for (int j = 0; j < len_pattern; j++){
			pattern += to_string(int(population.get_Mpop_geno()(i,j)));
		}
		list_sol.push_back(pattern);
	}


	vector<vector<string>> occ_pattern;
	int max_occ = 0;
	string best_pattern;
	int cpt;
	occ_pattern.resize(list_pattern.size(), vector<string>(2, ""));
	for (int i = 0; i < int(list_pattern.size()); i++){
		occ_pattern[i][0] = list_pattern[i];
		cpt = count(list_sol.begin(), list_sol.end(), list_pattern[i]);
		occ_pattern[i][1] = to_string(cpt);
		if (cpt > max_occ){
			max_occ = cpt;
			best_pattern = list_pattern[i];
		}
		else if (cpt == max_occ){
			best_pattern += (" and " + list_pattern[i]);
		}
		cout << "pattern " << occ_pattern[i][0] << " : " << occ_pattern[i][1] << endl;
	}
	cout << "The best pattern(s) is (are) " << best_pattern << " with a frequency of " << max_occ << endl;
	cout << "THE END 12";
}

//REFAIRE AVEC VECTOR DE VECTOR SUR LA FIN
//COMMENTER LE CODE GROS GUIGNOLE
