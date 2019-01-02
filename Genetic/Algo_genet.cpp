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
#include <string>
#include "Parametersfileparsing.hpp"
#include "Population.h"
#include "Parent.h"
#include "ContingencyTable.h"
#include "TheoricalTable.h"
#include "G2test.h"
#include "Child.h"
#include "Output.h"
#include "datainput.hpp"


using namespace std;


int main (int argc, char *argv[]) {
	srand (time(NULL));

	Parameters_file_parsing parameters("/home/courtin/Documents/M2/ProjetC/detection-of-epistasia-pattern/PARAMETERS_GENETIC.txt");


	char sep = parameters.sep;
	int nb_line_header = parameters.header_nrows;
	int len_pop = parameters.len_pop;
	int len_pattern = parameters.len_pattern;
	int nb_parents = parameters.nb_parents;
	float alpha = parameters.alpha;
	int nb_it = parameters.nb_it;
	float P_mutation = parameters.P_mutation;
	float P_selection = parameters.P_selection;
	string namefile = argv[1];
	string geno_path = argv[2];
	string pheno_path = argv[3];


	Data_input datas_geno(geno_path, sep, nb_line_header);
	int_matrix_type Mgeno = datas_geno.read();

	Data_input datas_pheno(pheno_path, sep, nb_line_header);
	int_matrix_type Mpheno = datas_pheno.read();

	Data_input header_line(geno_path, sep, nb_line_header);
	vector<string> header = header_line.get_snps();

	cout << Mgeno << endl;
	Population population(Mgeno, Mpheno, len_pop, len_pattern);

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


	float median;

	sort(G2_res.begin(), G2_res.end());

	median = G2_res[G2_res.size()/2]; //Median of the solutions' G2



	int iterator = 0;
	while (iterator < nb_it){

		Parent parents(len_pop, nb_parents, len_pattern, median, P_selection, population.get_Mpop_geno());
		parents.parents_selection();
		//parents.display_parents();

		Child children(population.get_Mpop_geno(), parents.get_MParents(), len_pattern, P_mutation, Mgeno.size2());
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
		cout << iterator << "/" << endl;
	}


	Output output(population.get_Mpop_geno(), header, len_pattern, len_pop, namefile);
	output.set_list_pattern();
	output.set_list_sol();
	output.set_best_sol();
	output.write_best_sol();

	cout << "THE END 12";
}
